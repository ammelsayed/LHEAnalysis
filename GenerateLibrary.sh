# 1) Generate dictionary
rootcling -f Particle_cxx.cxx -c Particle.h

# 2) Compile shared object
g++ -fPIC -shared -o Particle.so Particle.cxx Particle_cxx.cxx `root-config --cflags --libs`