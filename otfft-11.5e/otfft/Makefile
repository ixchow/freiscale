UNAME = ${shell uname}

ifeq ($(UNAME),Darwin)
    include Makefile.macOS
else
ifeq ($(UNAME),Linux)
    include Makefile.linux
else
ifeq ($(OS),Windows_NT)
    include Makefile.cygwin
else
    include Makefile.linux
endif
endif
endif
