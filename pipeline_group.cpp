//
// Created by markus on 01.10.19.
//

#include "pipeline_group.h"

void pipeline_group::control_pipeline() {
    // execute last pipeline stage for the next element
    pipeline_stage(stages - 1);

    // random sifting satisfied? automorphism probing satisfied?
    // -> done!

    // random sifting not satisfied yet?
    // -> generate random element and add to pipeline
    // automorphisms left?
    // -> add them to pipeline

    // keep track of how many elements are in pipeline, do not fill too much
}

void pipeline_group::pipeline_stage(int n) {
    // perform pipeline stage n
}
