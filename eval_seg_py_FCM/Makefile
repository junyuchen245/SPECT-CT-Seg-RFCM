CFLAGS=-g3 -O -fPIC -DDEBUG -std=c99

libsegment.so : miputil.o update_labels_new.o
	$(CC) $(CFLAGS) -shared -o libsegment.so miputil.o update_labels_new.o
