.DEFAULT_GOAL := help

GO ?= go
GOFMT ?= gofmt
PYTHON ?= python3

GO_DIRS := cmd codegen contracts topology
REFERENCE_DIR := references/benzene_toluene_5x2

.PHONY: help build test vet fmt fmt-check validate check test-debug reference-check

help:
	@printf '%s\n' \
	  'make build            Build CLI in bin/lumo' \
	  'make test             Run Go and generated C tests' \
	  'make vet              Run go vet' \
	  'make fmt              Format Go sources' \
	  'make fmt-check        Check Go formatting without changing files' \
	  'make validate         Validate both C99 builds against the 5x2 reference' \
	  'make check            Run fmt-check, vet, test and validate' \
	  'make test-debug       Run codegen tests and keep C debug artifacts' \
	  'make reference-check  Recompute and verify the independent reference'

build:
	mkdir -p bin
	$(GO) build -o bin/lumo ./cmd/lumo

test:
	$(GO) test -count=1 ./...

vet:
	$(GO) vet ./...

fmt:
	$(GOFMT) -w $(GO_DIRS)

fmt-check:
	@files=$$($(GOFMT) -l $(GO_DIRS)) || exit $$?; \
	if [ -n "$$files" ]; then \
	  printf 'Run make fmt for these files\n%s\n' "$$files"; \
	  exit 1; \
	fi

validate:
	$(PYTHON) $(REFERENCE_DIR)/validate.py --cc "$(CC)" --output build/validation

check: fmt-check vet test validate

test-debug:
	LUMO_CODEGEN_KEEP=1 $(GO) test -count=1 -v ./codegen/...

reference-check:
	$(PYTHON) $(REFERENCE_DIR)/calculate.py --check
