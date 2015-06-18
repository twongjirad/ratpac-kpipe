# -*- mode: Makefile -*-

#CC = clang++
CC = g++

ROOT_INCLUDES = `root-config --cflags`
ROOT_LIBS = `root-config --glibs` -lGeom
RAT_INCLUDES = -I./../include
RAT_LIBS = -L./../lib -lRATEvent -lrat
LOCAL_INC = -I. 


#CFLAGS = -g -fPIC 
CFLAGS = -O3 -fPIC
#CFLAGS += -std=c++11
CFLAGS += $(ROOT_INCLUDES) $(RAT_INCLUDES) $(LOCAL_INC)
LDFLAGS += `root-config --ldflags` 
LDLIBS = $(RAT_LIBS) $(ROOT_LIBS) -lstdc++

all: scrape_data libkptrigger.so analyze_data aux_data test_geotools veto_dark_noise_study

scrape_data.o: scrape_data.cc
	$(CC) $(CFLAGS) -c -o $@ $^

aux_data.o: aux_data.cc
	$(CC) $(CFLAGS) -c -o $@ $^

kptrigger.o: kptrigger.cc kptrigger.h
	$(CC) $(CFLAGS) -c -o $@ kptrigger.cc

geotools.o: geotools.cc geotools.h
	$(CC) $(CFLAGS) -c -o $@ geotools.cc

pmtinfo.o: pmtinfo.cc pmtinfo.hh
	$(CC) $(CFLAGS) -c -o $@ pmtinfo.cc

kpdaq.o: kpdaq.cc kpdaq.h
	$(CC) $(CFLAGS) -c -o $@ kpdaq.cc

tree.o: tree3.cc tree3.h
	$(CC) $(CFLAGS) -c -o $@ tree3.cc

gen_dark_noise.o: gen_dark_noise.cc gen_dark_noise.hh
	$(CC) $(CFLAGS) -c -o $@ gen_dark_noise.cc

prefitz.o: prefitz.cc prefitz.hh
	$(CC) $(CFLAGS) -c -o $@ prefitz.cc

libkptrigger.so: kptrigger.o pmtinfo.o gen_dark_noise.o prefitz.o kpdaq.o tree.o geotools.o
	$(CC) $(LDFLAGS) -shared -o $@ $^ $(LDLIBS)
	cp $@ ../lib/

scrape_data: scrape_data.o libkptrigger.so
	$(CC) $(LDFLAGS) -o $@ scrape_data.o $(LDLIBS) -lkptrigger

aux_data: aux_data.o libkptrigger.so
	$(CC) $(LDFLAGS) -o $@ aux_data.o $(LDLIBS) -lkptrigger

analyze_data.o: analyze_data.cc
	$(CC) $(CFLAGS) -c -o $@ $^

analyze_data: analyze_data.o libkptrigger.so
	$(CC) $(LDFLAGS) -o $@ analyze_data.o $(LDLIBS) -lkptrigger

test_geotools.o: test_geotools.cc
	$(CC) $(CFLAGS) -c -o $@ $^

test_geotools: test_geotools.o libkptrigger.so
	$(CC) $(LDFLAGS) -o $@ test_geotools.o $(LDLIBS) -lkptrigger

veto_dark_noise_study.o: veto_dark_noise_study.cc
	$(CC) $(CFLAGS) -c -o $@ $^

veto_dark_noise_study: veto_dark_noise_study.o libkptrigger.so
	$(CC) $(LDFLAGS) -o $@ veto_dark_noise_study.o $(LDLIBS) -lkptrigger

clean: 
	rm *.o scrape_data libkptrigger.so analyze_data aux_data test_geotools veto_dark_noise_study

