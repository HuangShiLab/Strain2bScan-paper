# Reproduction targets for the Strain2bScan multi-species benchmark.
# Usage:
#   export STRAIN2BSCAN_BIN=/path/to/strain2bscan      # or have `strain2bscan` on PATH
#   make genomes samples build scaling accuracy
# All work happens inside $(WORKDIR); large data stays out of git (see .gitignore).

WORKDIR ?= work
STRAIN2BSCAN_BIN ?= strain2bscan
SCRIPTS := $(CURDIR)/scripts
export WORKDIR STRAIN2BSCAN_BIN

$(WORKDIR):
	mkdir -p $(WORKDIR)

genomes: $(WORKDIR)              ## download 40 species x 4 strains from NCBI
	cd $(WORKDIR) && python3 $(SCRIPTS)/dl_multispecies.py

samples: $(WORKDIR)             ## simulate 20 multi-species samples (seeded) + truth
	cd $(WORKDIR) && python3 $(SCRIPTS)/sim_samples.py 20 8

build scaling:                  ## build per-species DBs + species/sample scaling gradients
	cd $(WORKDIR) && WORKDIR=$(WORKDIR) bash $(SCRIPTS)/run_ms.sh

accuracy:                       ## species/strain detection accuracy (default gate)
	cd $(WORKDIR) && python3 $(SCRIPTS)/ms_eval.py

gate-sweep:                     ## accuracy vs species-gate threshold
	cd $(WORKDIR) && for g in 30 100 200 400 800; do python3 $(SCRIPTS)/ms_eval2.py $$g; done

# ---- Reference-genome quality experiment (supplementary figure) ----
refqual-data: $(WORKDIR)        ## download C. acnes 64-panel + 5 mock reads + E. coli contaminant
	cd $(WORKDIR) && python3 $(SCRIPTS)/dl_acnes_panel.py

refqual: refqual-data           ## Strain2bScan arm: degrade DB across a quality gradient + figure
	cd $(WORKDIR) && python3 $(SCRIPTS)/run_refqual.py && python3 $(SCRIPTS)/plot_refqual.py
	@echo "outputs: $(WORKDIR)/refqual/refqual_degradation.tsv, $(WORKDIR)/figures/refqual_figure.*"

refqual-strainscan:             ## StrainScan arm (LINUX; needs $$STRAINSCAN, $$STRAINSCAN_PY)
	cd $(WORKDIR) && python3 $(SCRIPTS)/run_refqual_strainscan.py && python3 $(SCRIPTS)/plot_refqual.py

refqual-checkm2:                ## CheckM2 validation of simulated degraded genomes (LINUX)
	cd $(WORKDIR) && python3 $(SCRIPTS)/checkm2_validate.py

clean:                          ## remove the work directory
	rm -rf $(WORKDIR)

.PHONY: genomes samples build scaling accuracy gate-sweep \
        refqual-data refqual refqual-strainscan refqual-checkm2 clean
