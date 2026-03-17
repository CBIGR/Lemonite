#!/usr/bin/env bash

# LemonTree post-clustering analysis script
# Runs tight clustering and regulator analysis

# Parameters
REGULATOR_TYPES=$1
LEMONTREE_JAR=$2
WORK_DIR=$3
ORGANISM=$4         # e.g., "human", "mouse", "mmu"

cd ${WORK_DIR}

# Auto-replace TF list for mouse organism
if [[ "$ORGANISM" == "mouse" || "$ORGANISM" == "mmu" || "$ORGANISM" == "mus_musculus" ]]; then
    REGULATOR_TYPES=$(echo "$REGULATOR_TYPES" | sed -E 's/[Ll]overing_TF_list\.txt/lovering_TF_list_mouse.txt/g')
    echo "🐭 Organism is mouse - auto-replacing TF list with mouse version"
    echo "   Updated regulator_types: $REGULATOR_TYPES"
fi

# Check if required files exist
if [ ! -f "Preprocessing/LemonPreprocessed_expression.txt" ]; then
    echo "Error: Preprocessing/LemonPreprocessed_expression.txt not found"
    exit 1
fi

if [ ! -f "Preprocessing/LemonPreprocessed_complete.txt" ]; then
    echo "Error: Preprocessing/LemonPreprocessed_complete.txt not found"
    exit 1
fi

echo "====================================================="
echo "LemonTree Post-Clustering Analysis"
echo "====================================================="
echo "Regulator types: $REGULATOR_TYPES"
echo "LemonTree JAR: $LEMONTREE_JAR"
echo "Working directory: $WORK_DIR"
echo "====================================================="

# Build classpath including all dependencies
LEMONTREE_DIR=$(dirname ${LEMONTREE_JAR})
CLASSPATH="${LEMONTREE_JAR}"
if [ -d "${LEMONTREE_DIR}/lib" ]; then
    for jar in ${LEMONTREE_DIR}/lib/*.jar; do
        CLASSPATH="${CLASSPATH}:${jar}"
    done
fi

echo "Running LemonTree post-clustering analysis..."

# Create Lemon_out directory if it doesn't exist
mkdir -p Lemon_out

# First, we need to create a consensus clusterfile from all individual cluster outputs
# This step would typically be done by combining all cluster_* files into a single clusterfile
echo "Creating consensus cluster file..."
# Generate clusterfile with proper paths pointing to cluster files
rm -f Lemon_out/clusterfile
cluster_count=0

# Dynamically find all cluster files in Lemon_results directory
echo "Searching for cluster files..."
if [ -d "Lemon_results" ]; then
    for cluster_file in Lemon_results/cluster_*; do
        if [ -f "$cluster_file" ]; then
            echo "$cluster_file" >> Lemon_out/clusterfile
            cluster_count=$((cluster_count + 1))
        fi
    done
fi

# Also check in Lemon_out/Lemon_results if copied there
if [ -d "Lemon_out/Lemon_results" ]; then
    for cluster_file in Lemon_out/Lemon_results/cluster_*; do
        if [ -f "$cluster_file" ]; then
            # Convert path to relative from Lemon_results
            relative_path=$(echo "$cluster_file" | sed 's|Lemon_out/||')
            echo "$relative_path" >> Lemon_out/clusterfile
            cluster_count=$((cluster_count + 1))
        fi
    done
fi

echo "Found $cluster_count cluster files"

# If no cluster files found, create a basic clusterfile for testing
if [ ! -f "Lemon_out/clusterfile" ] || [ ! -s "Lemon_out/clusterfile" ]; then
    echo "No cluster files found, creating basic clusterfile for testing..."
    for i in {1..5}; do
        echo "Lemon_results/cluster_$i" >> Lemon_out/clusterfile
    done
fi

# Step 1: Generate tight clusters
echo "Generating tight clusters..."
java -cp ${CLASSPATH} lemontree.modulenetwork.RunCli \
    -task tight_clusters \
    -data_file Preprocessing/LemonPreprocessed_expression.txt \
    -cluster_file Lemon_out/clusterfile \
    -output_file Lemon_out/tight_clusters.txt \
    -node_clustering false \
    -min_weight 0.25 \
    -min_clust_size 10 \
    -min_clust_score 2

if [ $? -ne 0 ]; then
    echo "Error: tight_clusters task failed"
    exit 1
fi

# Step 2: Generate regulators for each regulator type dynamically
echo "Generating regulators for all regulator types..."

# Parse regulator_types parameter (format: "Prefix1:File1[:DataType1],Prefix2:File2[:DataType2],...")
# DataType is optional: 'c' for continuous (default), 'd' for discrete/binary
IFS=',' read -ra REGULATOR_PAIRS <<< "$REGULATOR_TYPES"

for pair in "${REGULATOR_PAIRS[@]}"; do
    PREFIX=$(echo "$pair" | cut -d':' -f1 | xargs)  # xargs trims whitespace
    FILENAME=$(echo "$pair" | cut -d':' -f2 | xargs)
    DATA_TYPE=$(echo "$pair" | cut -d':' -f3 | xargs 2>/dev/null)
    if [ -z "$DATA_TYPE" ]; then
        DATA_TYPE="c"
    fi
    
    # Convert prefix to lowercase for regulator list file and search multiple candidate locations
    PREFIX_LOWER=$(echo "$PREFIX" | tr '[:upper:]' '[:lower:]')
    
    # Determine if this is a TF list (contains "TF" in prefix or is a _list.txt file)
    # Discrete regulators are never TF lists — they contain abundance data
    IS_TF_LIST=false
    if [ "$DATA_TYPE" != "d" ]; then
        if [[ "$PREFIX" =~ [Tt][Ff] ]] || [[ "$FILENAME" =~ _list\.txt$ ]]; then
            IS_TF_LIST=true
        fi
    fi
    
    # TF lists: check data/ first (user override), then PKN/ (default)
    # Other regulators: only check Preprocessing/ and data/
    if [ "$IS_TF_LIST" = true ]; then
        CANDIDATE_FILES=( \
            "./Preprocessing/${PREFIX_LOWER}.txt" \
            "./LemonTree/Preprocessing/${PREFIX_LOWER}.txt" \
            "./data/${FILENAME}" \
            "./PKN/${FILENAME}" \
            "./${FILENAME}" \
        )
    else
        CANDIDATE_FILES=( \
            "./Preprocessing/${PREFIX_LOWER}.txt" \
            "./LemonTree/Preprocessing/${PREFIX_LOWER}.txt" \
            "./data/${FILENAME}" \
            "./${FILENAME}" \
        )
    fi

    REG_LIST_FILE=""
    for cand in "${CANDIDATE_FILES[@]}"; do
        if [ -f "$cand" ]; then
            REG_LIST_FILE="$cand"
            break
        fi
    done

    if [ -n "$REG_LIST_FILE" ]; then
        echo ""
        echo "Processing $PREFIX regulators..."
        echo "  Regulator list: $REG_LIST_FILE" 
        echo "  Data type: $DATA_TYPE (c=continuous, d=discrete)"
        echo "  Output prefix: Lemon_out/$PREFIX"

        java -cp ${CLASSPATH} lemontree.modulenetwork.RunCli \
            -task regulators \
            -data_file Preprocessing/LemonPreprocessed_complete.txt \
            -reg_file "$REG_LIST_FILE" \
            -cluster_file Lemon_out/tight_clusters.txt \
            -output_file Lemon_out/$PREFIX

        if [ $? -ne 0 ]; then
            echo "Error: $PREFIX regulators task failed"
            exit 1
        fi

        echo "✅ Successfully generated $PREFIX regulators (from $REG_LIST_FILE)"
    else
        echo "Warning: Regulator list not found in any candidate location for: ${PREFIX} (expected ${FILENAME}), skipping $PREFIX"
    fi
done

echo ""
echo "LemonTree post-clustering analysis completed successfully"
