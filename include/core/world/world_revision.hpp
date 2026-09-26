#pragma once

#include "util/typedefs.hpp"

using WorldRevision = u64;

struct WorldRevisions {
    WorldRevision topology = 0; // change body collections
    WorldRevision dynamics = 0; 
    WorldRevision providers = 0;
    WorldRevision instruments = 0;
    WorldRevision graphics = 0;
    WorldRevision state = 0;
};