
CXX       := g++
CXXFLAGS  := -O3 -std=c++17
TARGET    := tsp_sa
SRC       := tsp_sa.cpp
PLOT      := routeplot.py
SEED      := 1

DATA_150  := cities150.dat
DATA_1K   := cities1k.dat
DATA_2K   := cities2k.dat

.PHONY: all build run clean

# Default target
all: build run

build: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

run: run150 run1k run2k

run150: $(TARGET)
	./$(TARGET) $(DATA_150) $(SEED)
	python $(PLOT) $(DATA_150) best_route_150.dat \
		-s schedule_150.csv --save-png --no-show

run1k: $(TARGET)
	./$(TARGET) $(DATA_1K) $(SEED)
	python $(PLOT) $(DATA_1K) best_route_1k.dat \
		-s schedule_1k.csv --save-png --no-show

run2k: $(TARGET)
	./$(TARGET) $(DATA_2K) $(SEED)
	python $(PLOT) $(DATA_2K) best_route_2k.dat \
		-s schedule_2k.csv --save-png -w --no-show

clean:
	rm -f $(TARGET) \
	      best_route_*.dat \
	      schedule_*.csv \
	      an*.png \
	      *.o *.pdf
