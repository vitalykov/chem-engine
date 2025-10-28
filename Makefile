SRC := ./src
TEST_EXEC := $(SRC)/test

.PHONY: test
test: $(TEST_EXEC)
	cd $(SRC) && ./test

$(TEST_EXEC):
	g++ -Wall $(SRC)/*.cpp -o $@

.PHONY: clean
clean:
	rm -f $(TEST_EXEC)
