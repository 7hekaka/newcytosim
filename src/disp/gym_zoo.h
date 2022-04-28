// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

struct flute2;

namespace gym
{
    size_t zoo_init(flute2*, flute2* const);

    void zoo_stroke(char, float width);

    void zoo_paint(char);
};
