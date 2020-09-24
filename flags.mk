CXXFLAGS := -std=c++17 -m64 -pedantic -pedantic-errors -Werror -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual -Wformat=2 -Weffc++ -Wold-style-cast
CFLAGS := -m64 -pedantic -pedantic-errors -std=gnu99 -Werror -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual -Wformat=2 -Wstrict-prototypes -Wmissing-prototypes
DEBUG ?= 0
ifeq ($(DEBUG),1)
	DFLAGS := -g -Og -ggdb3
else
	DFLAGS := -O3 -flto
endif
CFLAGS += $(DFLAGS)
CXXFLAGS += $(DFLAGS)

# clang is dumb
ifneq ($(findstring clang,$(CXX)),)
	CXXFLAGS += -Wno-missing-braces -Wmissing-field-initializers
endif

override LINTFLAGS += -checks=performance-*,modernize-*,readability-*,misc-*,hicpp-*,cppcoreguidelines-*,google-*,-readability-magic-numbers,-readability-misleading-indentation,-hicpp-signed-bitwise,-hicpp-member-init,-cppcoreguidelines-avoid-magic-numbers,-cppcoreguidelines-pro-type-member-init,-cppcoreguidelines-pro-type-reinterpret-cast,-cppcoreguidelines-pro-bounds-pointer-arithmetic,-cppcoreguidelines-pro-bounds-constant-array-index,-google-runtime-references,-*braces-around-statements,-readability-avoid-const-params-in-decls -header-filter=.*
.PHONY: lintphony
lint_%: CXXFLAGS+= -Wno-missing-braces -Wmissing-field-initializers
lint_%: %.cc
	clang-tidy $< $(LINTFLAGS) -- $(CXXFLAGS) $(CPPFLAGS)
