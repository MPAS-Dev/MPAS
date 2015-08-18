.SUFFIXES: .c .o

EZXML_PATH= ../../external/ezxml

SHARED_OBJS = fortprintf.o
PARSE_OBJS = parse.o dictionary.o gen_inc.o utility.o
FILTER_OBJS = filter_parse.o

all: ezxml parse filter-parse
	($(MAKE) parse CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)")

ezxml:
	(cd $(EZXML_PATH); $(MAKE) CFLAGS="$(CFLAGS) $(TOOL_TARGET_ARCH)")

parse: ezxml $(PARSE_OBJS) $(SHARED_OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(EZXML_PATH)/ezxml.o  -I$(EZXML_PATH) -o $@ $(PARSE_OBJS) $(SHARED_OBJS)

filter-parse: ezxml $(FILTER_OBJS) $(SHARED_OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(EZXML_PATH)/ezxml.o -I$(EZXML_PATH) -o $@ $(FILTER_OBJS) $(SHARED_OBJS)

clean:
	$(RM) *.o parse
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.c.o: ezxml
	$(CC) $(TOOL_TARGET_ARCH) $(CFLAGS) $(CPPFLAGS) -I$(EZXML_PATH) -c $<
