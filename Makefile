TESTS_DIR = tests
TESTS_EXE = $(TESTS_DIR)/test-runner

PHONY: test
test:
	@cmake -B tests -S .
	@cmake --build $(TESTS_DIR)
	./$(TESTS_EXE)

.PHONY: clean
clean:
	rm -rf $(TESTS_DIR)
