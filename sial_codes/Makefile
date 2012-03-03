
SIAL_COMPILER=../bin/sial
SIO_DIR=../bin/sio

SIAL_SOURCE=$(wildcard *.sial)
SIO_OBJS=$(SIAL_SOURCE:.sial=.sio)

all: $(SIO_OBJS) $(SIO_DIR)/sial_config $(SIO_DIR)/default_jobflows

binclean clean distclean: % : ;
	rm -f *.sio

rebuild relink: % : binclean all

libclean ppclean archive: % : ;

$(SIO_DIR)/sial_config: sial_config
	cp sial_config $(SIO_DIR)

$(SIO_DIR)/default_jobflows: default_jobflows
	cp default_jobflows $(SIO_DIR)

.SUFFIXES: .sial .sio

.sial.sio:
	$(SIAL_COMPILER) $< && cp $@ $(SIO_DIR)
