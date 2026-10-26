# EMBD downstream analysis (Python 3 port): common tasks.
#
#   make setup              locked environment in .venv (needs uv)
#   make test               fast tests (~1 min)
#   make reproduce          spreadsheet -> every output in $(OUT) (~2.5 min, < 0.5 GB RAM)
#   make test-slow          compare every output in $(OUT) with tests/golden
#   make mem-profile        peak RSS and time per pipeline step -> $(OUT)_mem/mem_profile.tsv
#   make docker-build       runtime image (embd:latest) and test image (embd:test)
#   make docker-reproduce   full pipeline in the runtime image, memory capped at $(DOCKER_MEMORY)
#   make docker-test        fast + slow tests in the test image, against the outputs in $(OUT)
#   make apptainer          embd.sif from apptainer.def (HPC; needs apptainer)

OUT ?= out
DOCKER_MEMORY ?= 3g
PY := .venv/bin/python
SHEET := data/All_genes_mySQL_table.xlsx
INDEX := data/Index_EMBD_ID_Gene_names.xlsx
ABS_OUT = $(abspath $(OUT))

.PHONY: setup test reproduce test-slow mem-profile docker-build docker-reproduce docker-test apptainer clean

setup:
	uv sync --frozen --extra test
	-chflags nohidden .venv/lib/python3.12/site-packages/*.pth 2>/dev/null  # macOS, DECISIONS P2-8

test:
	$(PY) -m pytest -q

reproduce:
	$(PY) -m embd.pipeline reproduce --sheet $(SHEET) --index $(INDEX) --out $(OUT)

test-slow:
	EMBD_PORT_RUN_DIR=$(ABS_OUT) $(PY) -m pytest -q -m slow

mem-profile:
	tools/mem_profile.sh $(OUT)_mem

docker-build:
	docker build -t embd:latest .
	docker build --target test -t embd:test .

docker-reproduce:
	mkdir -p $(ABS_OUT)
	docker run --rm --memory=$(DOCKER_MEMORY) --memory-swap=$(DOCKER_MEMORY) --user $$(id -u):$$(id -g) \
		-v $(abspath data):/data:ro -v $(ABS_OUT):/out embd:latest

docker-test:
	docker run --rm --memory=$(DOCKER_MEMORY) --memory-swap=$(DOCKER_MEMORY) \
		-e EMBD_PORT_RUN_DIR=/run -v $(ABS_OUT):/run:ro embd:test -q -m "slow or not slow"

apptainer:
	apptainer build embd.sif apptainer.def

clean:
	rm -rf build $(OUT) $(OUT)_mem .pytest_cache
