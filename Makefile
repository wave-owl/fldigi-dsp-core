# fldigi-dsp-core -- plain make, C99, no dependencies besides libm.
#   make              static library into $(BUILD)
#   make test         build and run the golden tests
#   make sanitize     the tests under AddressSanitizer + UndefinedBehaviorSanitizer
#   make golden       regenerate tests/golden from fldigi 4.2.13 (FLDIGI=/path/to/fldigi)
#                     (Jalocha headers, and fldigi's own filter/FEC/varicode sources)
#   make checksum     print the content checksum that family cores pin in dsp.lock
#
# Family cores (fldigi-olivia-core, ...) consume this repository as the submodule dsp/
# and compile $(DSP_SOURCES) themselves after checking dsp.lock (R-920.3).
BUILD ?= build
CC ?= cc
CXX ?= c++
# -ffp-contract=off: no fused multiply-add, the arithmetic follows the reference order.
CFLAGS ?= -O3
CORE_FLAGS = -std=c99 -Wall -Wextra -Werror -pedantic -ffp-contract=off -fPIC -Iinclude
LDLIBS = -lm

DSP_SOURCES = src/fdsp_jalocha.c src/fdsp_gfft.c src/fdsp_filters.c src/fdsp_fec.c
DSP_HEADERS = include/fldigi_dsp.h
OBJ = $(patsubst src/%.c,$(BUILD)/%.o,$(DSP_SOURCES))

all: static
static: $(BUILD)/libfldigi_dsp.a

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.c $(DSP_HEADERS) | $(BUILD)
	$(CC) $(CORE_FLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/libfldigi_dsp.a: $(OBJ)
	ar rcs $@ $^

$(BUILD)/test_dsp: tests/test_dsp.c $(BUILD)/libfldigi_dsp.a
	$(CC) $(CORE_FLAGS) $(CFLAGS) $< $(BUILD)/libfldigi_dsp.a -o $@ $(LDLIBS)

test: $(BUILD)/test_dsp
	$(BUILD)/test_dsp tests/golden

SAN_FLAGS = -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all

sanitize:
	$(MAKE) --no-print-directory BUILD=$(BUILD)/san CFLAGS="$(SAN_FLAGS)" \
		LDLIBS="$(LDLIBS) $(SAN_FLAGS)" test

golden:
	@test -n "$(FLDIGI)" || { echo "set FLDIGI=/path/to/fldigi (tag v4.2.13)"; exit 2; }
	mkdir -p $(BUILD)
	$(CXX) -std=c++11 -O2 -ffp-contract=off -w -I$(FLDIGI)/src/include/jalocha \
		tests/golden/make_golden.cxx -o $(BUILD)/make_golden
	$(BUILD)/make_golden tests/golden
	$(CXX) -std=c++11 -O2 -ffp-contract=off -w -Itests/golden/stub -I$(FLDIGI)/src/include \
		tests/golden/make_golden_fldigi.cxx $(FLDIGI)/src/filters/fftfilt.cxx \
		$(FLDIGI)/src/filters/filters.cxx $(FLDIGI)/src/filters/viterbi.cxx \
		$(FLDIGI)/src/mfsk/interleave.cxx $(FLDIGI)/src/mfsk/mfskvaricode.cxx \
		$(FLDIGI)/src/misc/misc.cxx -o $(BUILD)/make_golden_fldigi
	$(BUILD)/make_golden_fldigi tests/golden

# Content checksum of the sources a family core links (sorted path list, then each
# file's bytes): the value recorded in a family core's dsp.lock.
checksum:
	@sh tools/dsp_checksum.sh .

clean:
	rm -rf $(BUILD)

.PHONY: all static test sanitize golden checksum clean
