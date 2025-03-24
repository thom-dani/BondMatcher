#!/bin/bash

# Dependencies

sudo apt install -y cmake ninja-build libboost-system-dev libopengl-dev libxcursor-dev
sudo apt install -y qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools qttools5-dev qtxmlpatterns5-dev-tools libqt5x11extras5-dev libqt5svg5-dev qttools5-dev libqt5x11extras5-dev libqt5svg5-dev qtxmlpatterns5-dev-tools
sudo apt install -y libcgal-dev
sudo apt install libsqlite3-dev
sudo apt install -y python3-sklearn python3-pip
sudo apt install -y git
sudo apt install build-essential

pip install umap-learn --break-system-packages

# ParaView installation

git clone https://github.com/topology-tool-kit/ttk-paraview.git
cd ttk-paraview
mkdir build
cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DCMAKE_INSTALL_PREFIX=../install ..
ninja install

PV_PREFIX=`pwd`/../install
export PATH=$PATH:$PV_PREFIX/bin
export LD_LIBRARY_PATH=$PV_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH:$PV_PREFIX/lib/python3.12/site-packages

# TTK installation

cd ../../ttk
mkdir build
cd build
paraviewPath=`pwd`/../../ttk-paraview/install/lib/cmake/paraview-5.13
cmake -G Ninja -DCMAKE_INSTALL_PREFIX=../install -DParaView_DIR=$paraviewPath -DTTK_ENABLE_SQLITE3=1 ..
ninja install

TTK_PREFIX=`pwd`/../install
export PV_PLUGIN_PATH=$TTK_PREFIX/bin/plugins/TopologyToolKit
export LD_LIBRARY_PATH=$TTK_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH:$TTK_PREFIX/lib/python3.12/site-packages

cd ../..
