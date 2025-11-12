TESTS_DIR = build
TESTS_EXE = $(TESTS_DIR)/test-runner

PHONY: test
test: $(TESTS_EXE)
	./$(TESTS_EXE)

$(TESTS_EXE):
	@cmake -B $(TESTS_DIR) -S .
	@cmake --build $(TESTS_DIR) --parallel

.PHONY: clean
clean:
	rm -rf $(TESTS_DIR)
