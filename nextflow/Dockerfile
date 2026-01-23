# LemonTree Analysis Pipeline Docker Image
FROM rocker/r-ver:4.3.0

# Set labels
LABEL maintainer="Boris Vandemoortele"
LABEL description="Docker image for Lemonite multi-omics analysis pipeline"
LABEL version="1.0.0"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONNOUSERSITE=1
ENV R_PROFILE_USER=/.Rprofile
ENV R_ENVIRON_USER=/.Renviron

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # System tools
    wget \
    curl \
    git \
    unzip \
    build-essential \
    # Java for LemonTree and Nextflow (requires Java 17+)
    openjdk-17-jdk \
    # Python
    python3 \
    python3-pip \
    python3-dev \
    # R compilation dependencies
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    # Additional system dependencies for R packages
    libgsl-dev \
    libmpfr-dev \
    libgmp-dev \
    libglpk-dev \
    # Bioconductor dependencies
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    # Network tools
    libigraph-dev \
    && rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64

# Install Python packages
RUN pip3 install --no-cache-dir \
    pandas==1.5.3 \
    numpy==1.24.3 \
    scipy==1.10.1 \
    matplotlib==3.7.1 \
    seaborn==0.12.2 \
    networkx==3.1 \
    scikit-learn==1.2.2 \
    statsmodels==0.14.0 \
    plotly==5.14.1 \
    megago \
    argparse

# Install R packages
RUN R -e "install.packages('BiocManager', repos='https://cran.rstudio.com')"

# Install CRAN packages (core packages first)
RUN R -e "install.packages(c( \
    'optparse', \
    'data.table', \
    'tidyverse', \
    'dplyr', \
    'ggplot2', \
    'pheatmap', \
    'stringr', \
    'parallel', \
    'ggrepel', \
    'enrichR', \
    'future.apply', \
    'future', \
    'doParallel', \
    'httr', \
    'jsonlite', \
    'plotly', \
    'htmlwidgets', \
    'DT', \
    'RColorBrewer', \
    'devtools', \
    'Hmisc', \
    'foreach' \
), repos='https://cran.rstudio.com', dependencies=TRUE, force=TRUE)"

# Install Bioconductor packages (including WGCNA dependencies)
RUN R -e "BiocManager::install(c( \
    'BiocGenerics', \
    'S4Vectors', \
    'IRanges', \
    'GenomeInfoDb', \
    'GenomicRanges', \
    'Biobase', \
    'SummarizedExperiment', \
    'DESeq2', \
    'biomaRt', \
    'decoupleR', \
    'OmnipathR', \
    'AnnotationDbi', \
    'fgsea', \
    'org.Hs.eg.db', \
    'org.Mm.eg.db', \
    'ReactomePA', \
    'BiocParallel', \
    'impute', \
    'preprocessCore', \
    'GO.db', \
    'KEGG.db', \
    'limma', \
    'edgeR', \
    'GSVA' \
), ask=FALSE, force=TRUE)"

# Install WGCNA separately (it has special dependencies)
RUN R -e "BiocManager::install('WGCNA', ask=FALSE, force=TRUE)" && \
    R -e "library(WGCNA); cat('WGCNA loaded successfully\n')"

# Download CollecTRI network during build to avoid runtime download issues on HPC
RUN mkdir -p /opt/lemonite_data && \
    R -e "library(decoupleR); \
    cat('Downloading CollecTRI network for offline use...\n'); \
    net <- tryCatch({ \
      get_collectri(organism='human', split_complexes=FALSE) \
    }, error = function(e) { \
      cat('Error downloading CollecTRI:', e\$message, '\n'); \
      return(NULL); \
    }); \
    if (!is.null(net)) { \
      saveRDS(net, '/opt/lemonite_data/collectri_network.rds'); \
      cat('CollecTRI network saved with', nrow(net), 'interactions\n'); \
    } else { \
      cat('WARNING: CollecTRI network download failed during build\n'); \
    }"

# Install igraph with system libraries
RUN R -e "install.packages('igraph', repos='https://cran.rstudio.com', dependencies=TRUE)"

RUN R -e "install.packages(c('tidytree', 'yulab.utils', 'aplot', 'ggfun', 'gson', 'ggplotify'), repos='https://cran.rstudio.com')"

# Install fgsea separately with proper dependencies
RUN R -e "BiocManager::install('fgsea', ask=FALSE, force=TRUE)"

# Install GOSemSim first (dependency for DOSE)
RUN R -e "BiocManager::install('GOSemSim', ask=FALSE, force=TRUE)"

# Install DOSE (depends on fgsea)
RUN R -e "BiocManager::install('DOSE', ask=FALSE, force=TRUE)"

# Install ggtree dependencies then ggtree
RUN R -e "BiocManager::install('treeio', ask=FALSE, force=TRUE)" && \
    R -e "BiocManager::install('ggtree', ask=FALSE, force=TRUE)"

# Install enrichplot (depends on DOSE and ggtree)
RUN R -e "BiocManager::install('enrichplot', ask=FALSE, force=TRUE)"

# Install ReactomePA
RUN R -e "BiocManager::install('ReactomePA', ask=FALSE, force=TRUE)"

# Install clusterProfiler
RUN R -e "BiocManager::install('clusterProfiler', ask=FALSE, force=TRUE)" && \
    R -e "library(clusterProfiler); cat('clusterProfiler loaded successfully\n')"

# Install remaining packages
RUN R -e "BiocManager::install(c('dorothea'), ask=FALSE, force=TRUE)"

# Create directories
RUN mkdir -p /opt/lemontree /app /data /results

# Copy LemonTree JAR file and dependencies from local system
COPY LemonTree/ /opt/lemontree/

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Copy pipeline files
COPY . /app/

# Copy PKN files
COPY PKN/ /opt/PKN/

# Set working directory
WORKDIR /app

# Set environment variables for LemonTree
ENV LEMONTREE_JAR=/opt/lemontree/lemontree_v3.1.1.jar

# Create entrypoint script
RUN echo '#!/bin/bash' > /app/entrypoint.sh && \
    echo 'set -e' >> /app/entrypoint.sh && \
    echo '' >> /app/entrypoint.sh && \
    echo '# Default command if no arguments provided' >> /app/entrypoint.sh && \
    echo 'if [ $# -eq 0 ]; then' >> /app/entrypoint.sh && \
    echo '    echo "LemonTree Analysis Pipeline Docker Container"' >> /app/entrypoint.sh && \
    echo '    echo "============================================="' >> /app/entrypoint.sh && \
    echo '    echo ""' >> /app/entrypoint.sh && \
    echo '    echo "Usage examples:"' >> /app/entrypoint.sh && \
    echo '    echo "  # Run pipeline:"' >> /app/entrypoint.sh && \
    echo '    echo "  docker run -v /path/to/data:/data lemontree-pipeline nextflow run main.nf --input_dir /data"' >> /app/entrypoint.sh && \
    echo '    echo ""' >> /app/entrypoint.sh && \
    echo '    echo "  # Interactive shell:"' >> /app/entrypoint.sh && \
    echo '    echo "  docker run -it -v /path/to/data:/data lemontree-pipeline bash"' >> /app/entrypoint.sh && \
    echo '    echo ""' >> /app/entrypoint.sh && \
    echo '    echo "Available commands: nextflow, Rscript, python3, java"' >> /app/entrypoint.sh && \
    echo '    exit 0' >> /app/entrypoint.sh && \
    echo 'fi' >> /app/entrypoint.sh && \
    echo '' >> /app/entrypoint.sh && \
    echo '# Execute the provided command' >> /app/entrypoint.sh && \
    echo 'exec "$@"' >> /app/entrypoint.sh && \
    chmod +x /app/entrypoint.sh

# Set the entrypoint
ENTRYPOINT ["/app/entrypoint.sh"]

# Default command
CMD ["bash"]
