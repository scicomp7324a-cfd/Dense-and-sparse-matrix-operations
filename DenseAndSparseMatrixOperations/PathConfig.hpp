# pragma once 

#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <map>
#include <functional>

using namespace std ;
namespace fs = filesystem ;

class PathConfig{

    private:
        fs::path baseDir_ ;
        fs::path outputDir_ ;

        using MemFn = fs::path (PathConfig::*)() const;
        map<int, MemFn> meshDirectoryMap ;
        map<int, MemFn> mtxDirectoryMap ;

    public:
        PathConfig()
            :baseDir_("Files"), outputDir_("Output"){
                createMeshDirectoryMap() ;
                createMtxDirectoryMap() ;
            }

        PathConfig(const fs::path& baseDir, const fs::path& outputDir)
            :baseDir_(baseDir), outputDir_(outputDir){
                createMeshDirectoryMap() ;
                createMtxDirectoryMap() ;
            }

        // Base and Output directories 
        const fs::path& baseDir() const {return baseDir_;}
        const fs::path& outputDir() const {return outputDir_;}

        // Sub folders
        fs::path meshDir() const{return baseDir_ / "MeshFiles" ; } 
        fs::path mtxDir() const{return baseDir_ / "MTXFiles" ; }
        
        // Dimension based sub folder 
        fs::path dim2x5Dir() const{return meshDir() / "2x5" ;}
        fs::path dim4x5Dir() const{return meshDir() / "4x5" ;}
        fs::path dim5x10Dir() const{return meshDir() / "5x10" ;}
        fs::path dim10x10Dir() const{return meshDir() / "10x10" ;}
        fs::path dim10x20Dir() const{return meshDir() / "10x20" ;}
        fs::path dim20x25Dir() const{return meshDir() / "20x25" ;}
        fs::path dim20x50Dir() const{return meshDir() / "20x50" ;}
        fs::path dim40x50Dir() const{return meshDir() / "40x50" ;}
        fs::path dim50x100Dir() const{return meshDir() / "50x100" ;}
        fs::path dim100x100Dir() const{return meshDir() / "100x100" ;}
        fs::path dim100x150Dir() const{return meshDir() / "100x150" ;}
        fs::path dim100x200Dir() const{return meshDir() / "100x200" ;}
        fs::path dim100x300Dir() const{return meshDir() / "100x300" ;}
        fs::path dim200x200Dir() const{return meshDir() / "200x200" ;}
        fs::path dim200x250Dir() const{return meshDir() / "200x250" ;}

        // Matrices sub foder 

        fs::path matrix10Dir() const{return mtxDir() / "matrices10" ;} 
        fs::path matrix20Dir() const{return mtxDir() / "matrices20" ;} 
        fs::path matrix50Dir() const{return mtxDir() / "matrices50" ;} 
        fs::path matrix100Dir() const{return mtxDir() / "matrices100" ;} 
        fs::path matrix200Dir() const{return mtxDir() / "matrices200" ;} 
        fs::path matrix500Dir() const{return mtxDir() / "matrices500" ;} 
        fs::path matrix1000Dir() const{return mtxDir() / "matrices1000" ;} 
        fs::path matrix2000Dir() const{return mtxDir() / "matrices2000" ;} 
        fs::path matrix5000Dir() const{return mtxDir() / "matrices5000" ;} 
        fs::path matrix10000Dir() const{return mtxDir() / "matrices10000" ;} 
        fs::path matrix15000Dir() const{return mtxDir() / "matrices15000" ;} 
        fs::path matrix20000Dir() const{return mtxDir() / "matrices20000" ;} 
        fs::path matrix30000Dir() const{return mtxDir() / "matrices30000" ;} 
        fs::path matrix40000Dir() const{return mtxDir() / "matrices40000" ;} 
        fs::path matrix50000Dir() const{return mtxDir() / "matrices50000" ;} 


        // Creation of a directory map
        void createMeshDirectoryMap(){
            meshDirectoryMap.emplace(10, &PathConfig::dim2x5Dir) ;
            meshDirectoryMap.emplace(20, &PathConfig::dim4x5Dir) ;
            meshDirectoryMap.emplace(50, &PathConfig::dim5x10Dir) ;
            meshDirectoryMap.emplace(100, &PathConfig::dim10x10Dir) ;
            meshDirectoryMap.emplace(200, &PathConfig::dim10x20Dir) ;
            meshDirectoryMap.emplace(500, &PathConfig::dim20x25Dir) ;
            meshDirectoryMap.emplace(1000, &PathConfig::dim20x50Dir) ;
            meshDirectoryMap.emplace(2000, &PathConfig::dim40x50Dir) ;
            meshDirectoryMap.emplace(5000, &PathConfig::dim50x100Dir) ;
            meshDirectoryMap.emplace(10000, &PathConfig::dim100x100Dir) ;
            meshDirectoryMap.emplace(15000, &PathConfig::dim100x150Dir) ;
            meshDirectoryMap.emplace(20000, &PathConfig::dim100x200Dir) ;
            meshDirectoryMap.emplace(30000, &PathConfig::dim100x300Dir) ;
            meshDirectoryMap.emplace(40000, &PathConfig::dim200x200Dir) ;
            meshDirectoryMap.emplace(50000, &PathConfig::dim200x250Dir) ;
        }

        void createMtxDirectoryMap(){
            mtxDirectoryMap.emplace(10, &PathConfig::matrix10Dir) ;
            mtxDirectoryMap.emplace(20, &PathConfig::matrix20Dir) ;
            mtxDirectoryMap.emplace(50, &PathConfig::matrix50Dir) ;
            mtxDirectoryMap.emplace(100, &PathConfig::matrix100Dir) ;
            mtxDirectoryMap.emplace(200, &PathConfig::matrix200Dir) ;
            mtxDirectoryMap.emplace(500, &PathConfig::matrix500Dir) ;
            mtxDirectoryMap.emplace(1000, &PathConfig::matrix1000Dir) ;
            mtxDirectoryMap.emplace(2000, &PathConfig::matrix2000Dir) ;
            mtxDirectoryMap.emplace(5000, &PathConfig::matrix5000Dir) ;
            mtxDirectoryMap.emplace(10000, &PathConfig::matrix10000Dir) ;
            mtxDirectoryMap.emplace(15000, &PathConfig::matrix15000Dir) ;
            mtxDirectoryMap.emplace(20000, &PathConfig::matrix20000Dir) ;
            mtxDirectoryMap.emplace(30000, &PathConfig::matrix30000Dir) ;
            mtxDirectoryMap.emplace(40000, &PathConfig::matrix40000Dir) ;
            mtxDirectoryMap.emplace(50000, &PathConfig::matrix50000Dir) ;
        }

        fs::path meshDirFor(int key) const {
            auto it = meshDirectoryMap.find(key);
            if (it == meshDirectoryMap.end()) {
                throw std::runtime_error("No mesh directory mapping for key = " + std::to_string(key));
            }
            return (this->*(it->second))(); // call the member function pointer
        }

        fs::path mtxDirFor(int key) const {
            auto it = mtxDirectoryMap.find(key);
            if (it == mtxDirectoryMap.end()) {
                throw std::runtime_error("No mtx directory mapping for key = " + std::to_string(key));
            }
            return (this->*(it->second))();
        }

        vector<string> getMeshFilenames(int dimension){

            fs::path meshPath = meshDirFor(dimension) ;
            vector<string> filenames = {meshPath / "points.txt",
                                        meshPath / "faces.txt",
                                        meshPath / "cells.txt",
                                        meshPath / "boundary.txt"} ;
            return filenames ;

        }

        vector<string> getMTXFilenames(int dimension){

            fs::path mtxPath = mtxDirFor(dimension) ;
            vector<string> filenames = {
                mtxPath / "laplacian_1.mtx",
                mtxPath / "laplacian_2.mtx",
                mtxPath / "laplacian_3.mtx",
                mtxPath / "laplacian_4.mtx"
            } ;

            return filenames ;
        }


        void validate() const {
            if (!fs::exists(baseDir_)) {
                throw std::runtime_error("Base directory not found: " + baseDir_.string());
            }
            if (!fs::exists(meshDir())) {
                throw std::runtime_error("Mesh directory not found: " + meshDir().string());
            }
            if (!fs::exists(mtxDir())) {
                throw std::runtime_error("Field directory not found: " + mtxDir().string());
            }
        }

        void ensureOutputDir() const {
            fs::create_directories(outputDir_);
        }

} ;

