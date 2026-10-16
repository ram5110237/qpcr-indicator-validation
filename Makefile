.PHONY: all data run test clean
all: data run test
data:
	Rscript scripts/simulate_example_data.R data/example
run:
	Rscript scripts/run_validation.R config/config.yaml
test:
	Rscript tests/testthat.R
clean:
	rm -rf results/example
