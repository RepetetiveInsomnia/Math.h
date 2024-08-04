CC= gcc
CFLAG= -Wall -Wextra -Werror -std=c11
CHECKFLAG= -lcheck
COVER=
COVER2=
r= 1
ifeq ($(shell uname),Linux)
		CHECKFLAG += -lsubunit -lrt -lm -lpthread
endif
ifeq ($(r), 1)
		COVER += -fprofile-arcs -ftest-coverage
		COVER2 += -fprofile-arcs
endif

all: s21_math.a

s21_math.o: s21_math.c
	$(CC) $(CFLAG) -c $^ -o $@

s21_math.a: s21_math.o
	ar rcs $@ $^
	rm -rf $^

s21_math_test.o: s21_math.c
	$(CC) $(CFLAG) $(COVER) -c $^ -o $@

s21_math_test.a: s21_math_test.o
	ar rcs $@ $^
	rm -rf $^

test_build: s21_math_test.a
	$(CC) $(CFLAG) $(COVER2) s21_test.c $^ $(CHECKFLAG) -o s21_test

test: test_build
	./s21_test || true

gcov_report: test
	lcov -c -d ./ -o s21_test.info
	genhtml -o ./report s21_test.info
	open ./report/index.html

clean:
	rm -rf *.o *.a s21_test.info s21_test *.gcda *gcno report