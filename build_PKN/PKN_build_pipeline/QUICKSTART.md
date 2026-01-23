# PKN Pipeline - Quick Start Guide

## Installation

### 1. Set Up Environment
```bash
cd /home/borisvdm/repo/LemonIte/build_PKN/python_scripts
bash setup.sh
```

This will:
- Create Python virtual environment
- Install all dependencies from requirements.txt
- Set up directory structure

### 2. Configure Paths
Edit `config.py` to point to your local database files:

```python
# Update these paths to match your system
HMDB_METABOLITES_XML = '/path/to/your/HMDB.xml'
BIOGRID_FILE = '/path/to/your/BIOGRID-ALL-4.4.238.tab3.txt'
STITCH_FILE = '/path/to/your/9606.protein_chemical.links.v5.0.tsv'
# ... etc
```

## Running the Pipeline

### Option 1: Run Complete Pipeline
```bash
python main.py --all
```
This runs all 3 steps sequentially (~3-5 hours total).

### Option 2: Run Individual Steps
```bash
# Step 1: Metabolite-gene interactions (~2-4 hours)
python main.py --step 1

# Step 2: PPI network (~30-60 minutes)
python main.py --step 2

# Step 3: Final integration (~5-10 minutes)
python main.py --step 3
```

### Option 3: Run Specific Databases Only
```bash
# Only query BioGRID, STITCH, and LINCS
python main.py --step 1 --databases biogrid,stitch,lincs
```

### Option 4: Resume from Checkpoint
```bash
# If L1000 or ChEMBL mapping was interrupted
python main.py --step 1 --resume
```

## Testing Individual Retrievers

Each retriever can be tested standalone:

```bash
cd step1_metabolites
python biogrid.py      # Test BioGRID retriever
python stitch.py       # Test STITCH retriever
python l1000.py        # Test L1000 retriever
# ... etc
```

## Expected Outputs

### Step 1 Output
- `output/metabolite_gene_PKN.tsv` - Combined metabolite-gene network
- `output/figures/database_upset_plot.png` - Database overlap visualization
- `output/figures/database_coverage_barplot.png` - Coverage per database

### Step 2 Output
- `output/PPI_network.tsv` - Combined protein-protein interaction network
- `output/figures/ppi_venn_diagram.png` - PPI source overlaps

### Step 3 Output
- `output/Final_PKN.tsv` - Unified metabolite-gene + PPI network
- `output/Final_PKN_with_links.tsv` - Network with database URLs
- `output/figures/database_comparison.png` - Coverage comparison charts
- `output/figures/network_statistics.png` - Degree distributions

## Monitoring Progress

### Check Logs
```bash
tail -f output/pipeline.log
```

### Check Output Files
```bash
ls -lh output/*.tsv
ls -lh output/figures/*.png
```

## Troubleshooting

### Problem: Missing Database Files
**Solution**: Update paths in `config.py` to point to correct locations

### Problem: API Rate Limits
**Solution**: Pipeline automatically retries with exponential backoff. Wait and it will resume.

### Problem: Out of Memory
**Solution**: 
- L1000 requires ~8GB RAM
- Close other applications
- Consider running on a machine with more memory

### Problem: L1000 Takes Too Long
**Solution**: 
- L1000 processes ~2.1GB file (1-2 hours is normal)
- Use `--resume` flag if interrupted
- Or exclude L1000: `python main.py --step 1 --databases biogrid,stitch,uniprot,...`

### Problem: ChEMBL Mapping Fails
**Solution**:
- Check internet connection (requires API access)
- Pipeline caches results to `HMDB_metabolites_ChEMBL_mapping.csv`
- Delete cache file to restart mapping

## Performance Tips

1. **Run overnight**: Full pipeline takes 3-5 hours
2. **Use --resume**: Checkpoints save progress for L1000 and ChEMBL mapping
3. **Check logs**: `pipeline.log` shows detailed progress
4. **Exclude slow databases**: Use `--databases` flag to skip L1000 if needed
5. **Run steps separately**: Split into 3 sessions if time is limited

## Directory Structure After Completion

```
python_scripts/
├── output/
│   ├── metabolite_gene_PKN.tsv           # Step 1 final
│   ├── PPI_network.tsv                   # Step 2 final
│   ├── Final_PKN.tsv                     # Step 3 final
│   ├── Final_PKN_with_links.tsv          # Step 3 with URLs
│   ├── BioGRID_interactions.csv          # Individual databases
│   ├── STITCH_interactions.csv
│   ├── ... (10 more database files)
│   ├── STRING_interactions.csv
│   ├── BioGRID_PPI_interactions.csv
│   ├── HuRI_interactions.csv
│   ├── HMDB_metabolites_ChEMBL_mapping.csv  # Cached mapping
│   ├── pipeline.log                      # Execution log
│   └── figures/
│       ├── database_upset_plot.png
│       ├── database_coverage_barplot.png
│       ├── ppi_venn_diagram.png
│       ├── database_comparison.png
│       └── network_statistics.png
```

## Next Steps

After pipeline completion:
1. **Review logs**: Check `pipeline.log` for any warnings
2. **Inspect outputs**: Open TSV files in Excel/pandas
3. **View visualizations**: Check PNG files in `output/figures/`
4. **Analyze network**: Use networkx or Cytoscape for further analysis

## Getting Help

- **Implementation details**: See `IMPLEMENTATION_SUMMARY.md`
- **Code documentation**: See docstrings in each module
- **Configuration**: See `config.py` for all settings
- **Examples**: See `__main__` blocks in each retriever file

## Citation

If you use this pipeline in your research, please cite:
- BioGRID: https://thebiogrid.org/
- STITCH: http://stitch.embl.de/
- STRING: https://string-db.org/
- UniProtKB: https://www.uniprot.org/
- IntAct: https://www.ebi.ac.uk/intact/
- chEMBL: https://www.ebi.ac.uk/chembl/
- LINCS: https://lincsproject.org/
- L1000: https://clue.io/
- Human-GEM: https://metabolicatlas.org/
- MetalinksDB: http://metalinks.csb.pitt.edu/
- HuRI: http://www.interactome-atlas.org/
