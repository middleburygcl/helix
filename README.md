### About

`helix` is a library that provides an interface for creating and manipulating hybrid-element meshes.

Currently, the main features include an interface for importing an `EGADS` geometry, optionally modifying it and then creating a volume mesh using `TetGen`. Additional meshing data structures will be added in the future.

#### Quickstart

Build with CMake:

```sh
git clone git@github.com:middpolymer/helix.git
cd helix
mkdir build
cd build
cmake ../
make
```

This will create a program called `helix` in the `bin` directory. There are currently two subprograms called `cad` and `mesh`. From your `build` directory, please run:

```sh
bin/helix cad --help
```

and 

```sh
bin/helix mesh --help
```

for information about how to use these programs. More detailed documentation for `helix` will be completed soon. 

### License

All `helix` source code (`.h`, `.cpp`) is distributed under the Apache-2.0 License.

Copyright 2023 Philip Claude Caplan

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


