
PROGRAM = bin/flIDE 

include make.sources

################################################################

OBJECTS = $(CFILES:.c=.o) $(CPPFILES:.cxx=.o) 

include makeinclude

.SUFFIXES : .cxx .c .o

.cxx.o :
	$(CXX) -Isrc -Iimages $(CXXFLAGS) -c $< -o $@
.c.o :
	$(CC) -Isrc -Iimages $(CFLAGS) -c $< -o $@

$(PROGRAM) : $(OBJECTS) 
	$(CXX) $(LDFLAGS) -o $(PROGRAM) $(OBJECTS) $(LDLIBS)

clean :
	-@ rm -f src/*.o $(PROGRAM) $(CLEAN) core *~ makedepend

depend:
	$(MAKEDEPEND) -Isrc -Iimages $(CXXFLAGS) $(CPPFILES) $(CFILES) > makedepend

makedepend:
	$(MAKEDEPEND) -Isrc -Iimages $(CXXFLAGS) $(CPPFILES) $(CFILES) > makedepend

include makedepend

install: $(PROGRAM)
	@echo "Installing FLUID..."
	@strip $(PROGRAM)
	@-mkdir -p $(bindir)
	@cp $(PROGRAM) $(bindir)/$(PROGRAM)
	@chmod 755 $(bindir)/$(PROGRAM)

uninstall:
	-@ rm -f $(bindir)/$(PROGRAM)

#
# Note: The rebuild target can only be used if you have the original .fl
#       files.  This is normally only used by the FLTK maintainers...
#

rebuild:
	./fluid -c about_panel.fl
	./fluid -c alignment_panel.fl
	./fluid -c function_panel.fl
	./fluid -c widget_panel.fl

#
# End of "$Id$".
#
