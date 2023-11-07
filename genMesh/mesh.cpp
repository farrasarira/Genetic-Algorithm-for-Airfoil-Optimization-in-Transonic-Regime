/*----------------------------------------------------------------------------------------------------------------------

        2D Airfoil mesher for Gmsh.

        Usage:

                mesh.exe bl_thickness hwall_n ratio hfar airfoil_list

        The first line is reserved for the name. The curve generated is a
        closed B-spline.  Points should be in a standard Selig format airfoil
        coordinate file.

        Usage example:

                mesh.exe 0.05 0.000001 1.3 1.5 "CH10-(smoothed).dat" e423.dat FX-84-W-150.dat S1223.dat


----------------------------------------------------------------------------------------------------------------------*/

#define _WIN32_WINNT 0x0500 // For console resizing

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
// #include <dos.h> //for delay

using namespace std;

int main(int argc, char *argv[])
{
  /* Resize Console to show presentation */


  /* Presentation */
  if (!argv[1] || argc < 5) {
    // MoveWindow(window_handle, x, y, width, height, redraw_window);
    cout << "\n------------------------------------------ 2D Airfoil mesher for Gmsh -------------------------------------------\n" << endl;
    cout << "\tUsage: \n\n\t\tmesh.exe bl_thickness hwall_n ratio hfar airfoil_list\n" << endl;
    cout << "\tThe first line is reserved for the name. The curve generated is a" << endl;
    cout << "\tclosed B-spline.  Points should be in a standard Selig format airfoil" << endl;
    cout << "\tcoordinate file." << endl;
    cout << "\n\tUsage example: \n\n\t\tmesh.exe 0.05 0.000001 1.3 1.5 \"CH10-(smoothed).dat\" e423.dat FX-84-W-150.dat S1223.dat\n" << endl;
    cout << "\n-----------------------------------------------------------------------------------------------------------------" << endl;
    return 0;
  } else {
    // MoveWindow(window_handle, x, y, width, height, redraw_window);
    cout << "\n------------------------------------------ 2D Airfoil mesher for Gmsh -------------------------------------------\n" << endl;
    cout << "\tMeshing: \n\n\t\t";
    for (int arg_i=5; arg_i<argc; arg_i++)
      cout << argv[arg_i] << ", ";
    cout << endl;
    cout << "\n\tMaximal thickness of the boundary layer: \n\n\t\t" << argv[1] << endl;
    cout << "\n\tMesh Size Normal to the The Wall: \n\n\t\t" << argv[2] << endl;
    cout << "\n\tSize Ratio Between Two Successive Layers: \n\n\t\t" << argv[3] << endl;
    cout << "\n\tElement size far from the wall: \n\n\t\t" << argv[4] << endl;
    cout << "\n\tGenerator function: hwall * ratio^(dist/hwall)\n" << endl;
    cout << "\n-----------------------------------------------------------------------------------------------------------------" << endl;
  }

  /* Boundary Layer Configuration Variables */
  float bl_thickness = atof(argv[1]);
  float bl_hwall_n = atof(argv[2]);
  float bl_ratio = atof(argv[3]);
  float bl_hfar = atof(argv[4]);

  /* Mesh airfoils in the list */
   for (int arg_i=5; arg_i<argc; arg_i++) {
     /* Read .dat airfoil file */
     vector<double> x,y;
     string x_coord, y_coord, airfoil_name;
     bool space1_flag, space2_flag, x_flag, y_flag, is_coord;
     ifstream airfoil_dat(argv[arg_i]);
     ofstream airfoil_geo;

     int curve_i = 0;
     int point_i = 0;
     // While the string continues
     for (string line; getline(airfoil_dat, line); ) {
       x_coord = y_coord = "";
       space1_flag = space2_flag = x_flag = y_flag = 0;
       if (point_i) { // Skip first line
         // For every position in the current line
         for (string::iterator it = line.begin(); it != line.end(); ++it) {
           string str(1,*it);
           is_coord = (str != " " && str != "\n");

           // Determine position in the line
           if ( !(space2_flag) && !(x_flag) && !(y_flag) && is_coord )
           x_flag = 1;
           if ( !(space2_flag) && x_flag && !(y_flag) && !(is_coord) )
           space2_flag = 1;
           if ( space2_flag && x_flag && !(y_flag) && is_coord )
           y_flag = 1;

           // Register coordinates
           if ( !(space2_flag) && x_flag && !(y_flag) && is_coord )
           x_coord = x_coord + *it;
           if ( space2_flag && x_flag && y_flag && is_coord )
           y_coord = y_coord + *it;
         }
         if (x_flag && y_flag)
         x.push_back(stod(x_coord));
         y.push_back(stod(y_coord));
       } else {
         airfoil_name = line;
         airfoil_geo.open(airfoil_name + ".geo");
         airfoil_geo << "Geometry.LineNumbers = 1; Geometry.SurfaceNumbers = 1;" << endl;
         //airfoil_geo << "Geometry.PointNumbers = 1;" << endl;
       }
       if (point_i == 0)
       point_i++;
     }
     x.pop_back();y.pop_back(); // erase (0,0) at the end

     /* Airfoil Points Creation in Gmsh File */
     for (int i =  0; i <= x.size(); i++)
      airfoil_geo << "Point(" << point_i++ << ") = {" << x[i] << "," << y[i] << ",0};" << endl;
     //cout << x[i]  << "," << y[i] << endl;

     // Points corresponding to the different sections
     int sp_start = 1;
     int sp_end = point_i-1;
     int le_up = floor(0.5*0.75*(sp_end-sp_start)+sp_start);
     int sp_50 = floor(0.5*(sp_end-sp_start)+sp_start);
     int le_do = floor(0.625*(sp_end-sp_start)+sp_start);

     /* Airfoil Curve Creation */
     int curve_loop_i = 0;
     if (x[x.size()] == x[0] && y[y.size()] == y[0]) { // If closed trailing edge
       airfoil_geo << "BSpline(" << ++curve_i << ") = {" << sp_start << ":" << le_up << "};" << endl;
       airfoil_geo << "Transfinite Curve {-" << curve_i << "} = 150 Using Progression 1.01;" << endl;
       airfoil_geo << "BSpline(" << ++curve_i << ") = {" << le_up << ":" << sp_50 << "};" << endl;
       airfoil_geo << "Transfinite Curve {-" << curve_i << "} = 50 Using Bump 1.01;" << endl;
       airfoil_geo << "BSpline(" << ++curve_i << ") = {" << sp_50 << ":" << le_do << "};" << endl;
       airfoil_geo << "Transfinite Curve {" << curve_i << "} = 50 Using Bump 1.01;" << endl;
       airfoil_geo << "BSpline(" << ++curve_i << ") = {" << le_do << ":" << sp_end << ",1};" << endl;
       airfoil_geo << "Transfinite Curve {" << curve_i << "} = 150 Using Progression 1.01;" << endl;
       airfoil_geo << "Curve Loop(" << ++curve_loop_i << ") = {" << curve_i-3 << "," << curve_i-2 << "," << curve_i-1 << "," << curve_i << "};" << endl;
     } else {
       airfoil_geo << "BSpline(" << ++curve_i << ") = {" << sp_start << ":" << le_up << "};" << endl;
       airfoil_geo << "Transfinite Curve {-" << curve_i << "} = 150 Using Progression 1.01;" << endl;
       airfoil_geo << "BSpline(" << ++curve_i << ") = {" << le_up << ":" << sp_50 << "};" << endl;
       airfoil_geo << "Transfinite Curve {-" << curve_i << "} = 50 Using Bump 1.01;" << endl;
       airfoil_geo << "BSpline(" << ++curve_i << ") = {" << sp_50 << ":" << le_do << "};" << endl;
       airfoil_geo << "Transfinite Curve {" << curve_i << "} = 50 Using Bump 1.01;" << endl;
       airfoil_geo << "BSpline(" << ++curve_i << ") = {" << le_do << ":" << sp_end << "};" << endl;
       airfoil_geo << "Transfinite Curve {" << curve_i << "} = 150 Using Progression 1.01;" << endl;
       airfoil_geo << "Line(" << ++curve_i << ") = {" << sp_end << "," << sp_start << "};" << endl; // Close Trailing Edge
       airfoil_geo << "Curve Loop(" << ++curve_loop_i << ") = {" << curve_i-4 << "," << curve_i-3 << "," << curve_i-2 << "," << curve_i-1 << "," << curve_i << "};" << endl;
     }

     /* Control Volume */
     airfoil_geo << "Point(" << point_i +1 << ") = {0.5,100,0};" << endl;
     airfoil_geo << "Point(" << point_i +2 << ") = {0.5,0,0};" << endl;
     airfoil_geo << "Point(" << point_i +3 << ") = {0.5,-100,0};" << endl;
     airfoil_geo << "Circle(" << ++curve_i << ") = {" << point_i +1 << "," << point_i +2 << "," << point_i +3 << "};" << endl;
     airfoil_geo << "Circle(" << ++curve_i << ") = {" << point_i +3 << "," << point_i +2 << "," << point_i +1 << "};" << endl;
     airfoil_geo << "Curve Loop(" << ++curve_loop_i << ") = {" << curve_i-1 << "," << curve_i << "};" << endl;

     /* Surface Creation */
     int surface_i = 0;
     airfoil_geo << "Plane Surface(" << ++surface_i << ") = {" << curve_loop_i << "," << curve_loop_i-1 << "};" << endl;

     /* Boundary Layer */
     airfoil_geo << "Field[1] = BoundaryLayer;" << endl;
     airfoil_geo << "Field[1].thickness = " << bl_thickness << ";" << endl;
     airfoil_geo << "Field[1].hwall_n = " << bl_hwall_n << ";" << endl;
     airfoil_geo << "Field[1].ratio = " << bl_ratio << ";" << endl;
     airfoil_geo << "Field[1].hfar = " << bl_hfar << ";" << endl;
     airfoil_geo << "Field[1].Quads = 1;" << endl;
     if (x[x.size()] == x[0] && y[y.size()] == y[0]) { // If closed trailing edge
       airfoil_geo << "Field[1].EdgesList = {1,2,3,4};" << endl;
       airfoil_geo << "Field[1].FanNodesList = {1};" << endl;
     } else {
       airfoil_geo << "Field[1].EdgesList = {1,2,3,4,5};" << endl;
       airfoil_geo << "Field[1].FanNodesList = {1," << x.size()+1 << "};" << endl;
     }
     airfoil_geo << "BoundaryLayer Field = 1;" << endl;

     /* Boundary Conditions */
     if (x[x.size()] == x[0] && y[y.size()] == y[0]) { // If closed trailing edge
       airfoil_geo << "Physical Curve(\"airfoil\") = {1, 2, 3, 4};" << endl;
       airfoil_geo << "Physical Curve(\"farfield\") = {5, 6};" << endl;
       airfoil_geo << "Physical Surface(\"fluid\") = {1};" << endl;
     } else {
       airfoil_geo << "Physical Curve(\"airfoil\") = {1, 2, 3, 4, 5};" << endl;
       airfoil_geo << "Physical Curve(\"farfield\") = {6, 7};" << endl;
       airfoil_geo << "Physical Surface(\"fluid\") = {1};" << endl;
     }

     /* Mesh */
     airfoil_geo << "SetOrder 2;" << endl; // Improves agreement of the bspline and the mesh
	   airfoil_geo << "Mesh.Format = 42;" << endl; // SU2 Format
     airfoil_geo << "Mesh 2;" << endl;
     airfoil_geo << "Save \"" << airfoil_name << ".su2\";" << endl; // Save 1D Mesh

     airfoil_dat.close();
     airfoil_geo.close();
    //  system(("gmsh \"" + airfoil_name + ".geo\"").c_str()); // Execute Gmsh with input script
    //  system(exit);
    //  slepp(1500);
    //  std::cout << "hallo" ;
   }
   return 0;
}

