CURRENT=`pwd`
CC	=g++
CFLAGS	=-O3 -g3 -pthread -Wall -c -fmessage-length=0
OTHERS	=-DHAVE_CONFIG_H 
EXP	=-rdynamic
OBJECTS	=$(CURRENT)/MIRZA.o $(CURRENT)/Parameters.o $(CURRENT)/Sequence.o  $(CURRENT)/Hybridize.o


all: MIRZA.o Parameters.o Sequence.o Hybridize.o
	$(CC) $(EXP) -o MIRZA $(OBJECTS)
Hybridize.o: ./src/Hybridize.cpp
	$(CC) $(CFLAGS) $(OTHERS) $(CURRENT)/src/Hybridize.cpp
Parameters.o: ./src/Parameters.cpp
	$(CC) $(CFLAGS) $(OTHERS) $(CURRENT)/src/Parameters.cpp
Sequence.o: ./src/Sequence.cpp
	$(CC) $(CFLAGS) $(OTHERS) $(CURRENT)/src/Sequence.cpp
MIRZA.o: ./src/MIRZA.o
	$(CC) $(CFLAGS) $(OTHERS) $(CURRENT)/src/MIRZA.cpp
clean:
	rm -rf *.o
	rm -rf MIRZA



