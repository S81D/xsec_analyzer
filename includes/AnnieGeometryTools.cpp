/// Set of tool for plotting working working ANNIE's geometry
// Aurther :: christian nguyen 

#include "AnnieGeometryTools.hh"



 // MicroBooNE
 double FV_X_MIN =   21.5;
 double FV_X_MAX =  234.85;

 double FV_Y_MIN = -95.0;
 double FV_Y_MAX =  95.0;

 double FV_Z_MIN =   21.5;
 double FV_Z_MAX =  966.8;
//
////ANNIE FV Boundaries
 double FV_RAD = 100.0;
//
 double A_FV_Y_MIN = -150.0;
 double A_FV_Y_MAX =  150.0;
//
 double A_Z_CTR = 168.1; //Z center of tank




TGraph* createBoundaryGraph(const std::string& plotType) {
    double radius = 100.;  // cm
    double y_min = -150.;  // cm
    double y_max = 150.;   // cm
    double z_center = 168.1;  // cm

    std::vector<Vertex_XYZ> VertexPostion;

    int n_points = 10000;  // Number of points to plot

    // Generate points along the boundary
    for (int i = 0; i < n_points; ++i) {
        double angle = 2 * M_PI * i / n_points;
        double x = radius * std::cos(angle);
        double z = radius * std::sin(angle) + z_center;
        double y = y_min + (y_max - y_min) * i / n_points;

          if (y > y_min && y < y_max && radius > std::sqrt((z - z_center) * (z - z_center) + x * x) && (z - z_center < radius)) {
              Vertex_XYZ input{x,y,z};
              VertexPostion.push_back(input);
              std::cout<<"(x,y,z) = "<< "(" << x<< ", "<< y << ", "<< z<< " ) "<<std::endl;
              
            }
         }

        if (plotType == "X vs Y") {
            TGraph  *TGraph_output = Make_X_vs_Y_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "X vs Z") {
            TGraph  *TGraph_output = Make_X_vs_Z_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "Y vs Z") {
            TGraph  *TGraph_output = Make_Y_vs_Z_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "R vs Y") {
            TGraph  *TGraph_output = Make_R_vs_Y_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "R vs X") {
            TGraph  *TGraph_output = Make_R_vs_X_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "R vs Z") {
            TGraph  *TGraph_output = Make_R_vs_Z_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "RR vs Z") {
            TGraph  *TGraph_output = Make_RR_vs_Z_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
        else if (plotType == "RR vs X") {
            TGraph  *TGraph_output = Make_RR_vs_X_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
         else if (plotType == "RR vs Y") {
            TGraph  *TGraph_output = Make_RR_vs_Y_Tgraph_fromVector(VertexPostion);
            return TGraph_output;
        }
       else{TGraph  *TGraph_output =nullptr;\
       return TGraph_output;
       }
    
    
}/////end of function 
