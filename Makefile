SRC_DIR := ./src
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
TEST_EXEC := $(SRC_DIR)/test

.PHONY: test
test: $(TEST_EXEC)
	cd $(SRC_DIR) && ./test

$(TEST_EXEC): $(SOURCES)
	g++ -Wall $^ -o $@

.PHONY: clean
clean:
	rm -f $(TEST_EXEC)
