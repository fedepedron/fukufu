///////////////////////////////////////////////////////////////////////////////
/* Usage:                                                                     /
/       ./fukui.o sourcefile.log                                              /
/ The file.log should file a Gaussian (preferably G09) single-point           /
/ calculation with "# pop=full iop(3/33=4)" options.                          /
///////////////////////////////////////////////////////////////////////////////
/ Compilation command:                                                        /
/      g++ src_fukui.cpp -o fukui.o -std=c++11                                /
/ The option -std=c++11 is necesary for 'stof'.                              */
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>


// Reads the number of MOs, alpha electrons, and beta electrons.
// file_name: the name of the input file.   (in)
// n_orb    : the number of MOs.            (out)
// elec_a   : the number of alpha electrons.(out)
// elec_b   : the number of beta electrons. (out)
void read_n_orbs(std::string file_name, int *n_orb, int *elec_a, int *elec_b)
{
    int           word_pos; // Reads word position in line.
    std::string   temp_str, line_read;
    std::ifstream input_file;

    input_file.open(file_name);
    if (input_file.is_open())
    {
        for(line_read; std::getline(input_file, line_read);){
            word_pos = line_read.find("primitive gaussians");
            if (word_pos > 0) {
                // Reads the number of orbital functions.
                temp_str = line_read.substr(1, 6);
                sscanf(temp_str.c_str(), "%i", n_orb);

                // Reads the number of alpha and beta electrons.
                std::getline(input_file, line_read);
                temp_str = line_read.substr(1, 6);
                sscanf(temp_str.c_str(), "%i", elec_a);
                temp_str = line_read.substr(26, 6);
                sscanf(temp_str.c_str(), "%i", elec_b);
            };
        };
    } else {
        std::cout << "ERROR: Cannot open input file. Please check filename. \n";
    };
    input_file.close();
    return;
};
//End: read_n_orbs


// Reads the overlap matrix, the MO-coefficients and the energies.
// file_name : the name of the input file.            (in)
// n_functs  : the number of MOs.                     (in)
// s_matrix  : the overlap matrix.                    (out)
// ener_a    : a vector containing alpha MO energies. (out)
// ener_b    : a vector containing beta MO energies.  (out)
// c_matrix_a: the alpha MO coefficient matrix.       (out)
// c_matrix_b: the beta MO coefficient matrix.        (out)
void read_matrix(std::string file_name, int n_functs, double *s_matrix,
	             int *at_of_orb, double *ener_a, double *ener_b, 
	             double *c_matrix_a, double *c_matrix_b)
{
    double        temp_num   ;
                  // Indexes for matrixes.
    int           index_x, index_y, index_mat_a, index_mat_b,
                  // Blocks, columns and lines to read in the file.
                  block, n_blocks, column, n_columns, line, n_lines,
                  word_pos, temp_int, current_atom;
    std::string   line_read, temp_str;
    std::ifstream input_file ;

    word_pos   = 0;
    n_blocks = (n_functs - n_functs % 5) / 5; 
    input_file.open(file_name);

    if (input_file.is_open())
    {
        for(line_read; std::getline(input_file, line_read);)
        {  
            // Reads Overlap matrix
            word_pos = line_read.find("*** Overlap ***");
            if (word_pos > 0) 
            {
                for (block = 1; block < n_blocks + 1 ; block = block+1)
                {
                    // Reads a line with non-relevant information.
                    std::getline(input_file, line_read);
                    n_lines = n_functs - 5*(block-1);                   
                    for (line = 1 ; line < n_lines + 1 ; line = line+1)
                    {
                        std::getline(input_file, line_read);
                        std::replace(line_read.begin(), line_read.end(),
                                     'D', 'E');
                        n_columns = line;
                        if (line > 5) { n_columns = 5; };
                        for (column = 1; column < n_columns+1; column=column+1)
                        {
                            temp_str = line_read.substr(8 + (column-1)*14, 13);
                            sscanf(temp_str.c_str(), "%lf", &temp_num);

                            index_x     = line   + 5*(block-1);
                            index_y     = column + 5*(block-1);
                            index_mat_a = index_x + n_functs*(index_y - 1);
                            index_mat_b = index_y + n_functs*(index_x - 1);
                            s_matrix[index_mat_a-1] = temp_num;
                            s_matrix[index_mat_b-1] = temp_num;
                        };
                    };               
                }; 
            };
            // Ends overlap matrix.

            // Reads Alpha Coefficients matrix and energies.
            word_pos = line_read.find("Alpha Molecular Orbital Coefficients");
            if (word_pos > 0) 
            {
                for (block = 1; block < n_blocks + 1 ; block += 1)
                {   // Reads two lines with non-relevant information.
                    std::getline(input_file, line_read);
                    std::getline(input_file, line_read);

                    // Reads the energies.
                    std::getline(input_file, line_read);
                    for (column = 1; column < 6; column += 1)
                    {
                        temp_str = line_read.substr(22 + (column-1)*10, 9);
                        sscanf(temp_str.c_str(), "%lf", &temp_num);
                        
                        index_x     = column + 5*(block-1);
                        ener_a[index_x-1] = temp_num;
                    };


                    // Reads the orbital coefficients.
                    for (line = 1 ; line < n_functs + 1 ; line += 1)
                    {
                        std::getline(input_file, line_read);
                     
                        // Gets the atom identity for each basis function.
                        if (block < 2)
                        {
                            temp_str = line_read.substr(5, 4);
                            sscanf(temp_str.c_str(), "%i", &temp_int);
                            if (temp_int > 0){ current_atom = temp_int; };
                            at_of_orb[line-1] = current_atom;
                        };

                        //Gets the MO coefficients.
                        for (column = 1; column < 6; column += 1)
                        {
                            temp_str = line_read.substr(22 + (column-1)*10, 9);
                            sscanf(temp_str.c_str(), "%lf", &temp_num);

                            index_x     = line;
                            index_y     = column + 5*(block-1);
                            index_mat_a = index_x + n_functs*(index_y - 1);
                            c_matrix_a[index_mat_a-1] = temp_num;
                        };
                    };               
                }; 
            };
            // Ends alpha coefficients.
            
            // Reads Beta Coefficients matrix and energies.
            word_pos = line_read.find("Beta Molecular Orbital Coefficients");
            if (word_pos > 0) 
            {
                for (block = 1; block < n_blocks + 1 ; block += 1)
                {   // Reads two lines with non-relevant information.
                    std::getline(input_file, line_read);
                    std::getline(input_file, line_read);

                    // Reads the energies.
                    std::getline(input_file, line_read);
                    for (column = 1; column < 6; column += 1)
                    {
                        temp_str = line_read.substr(22 + (column-1)*10, 9);
                        sscanf(temp_str.c_str(), "%lf", &temp_num);

                        index_x     = column + 5*(block-1);
                        ener_b[index_x-1] = temp_num;
                    };

                    // Reads the orbital coefficients.
                    for (line = 1 ; line < n_functs + 1 ; line += 1)
                    {
                        std::getline(input_file, line_read);
                        for (column = 1; column < 6; column += 1)
                        {
                            temp_str = line_read.substr(22 + (column-1)*10, 9);
                            sscanf(temp_str.c_str(), "%lf", &temp_num);

                            index_x     = line;
                            index_y     = column + 5*(block-1);
                            index_mat_a = index_x + n_functs*(index_y - 1);
                            c_matrix_b[index_mat_a-1] = temp_num;
                        };
                    };               
                }; 
            };
            // Ends Beta coefficients.
        };
    } else {
        std::cout << "ERROR: Cannot open input file. Please check filename. \n";
    };
    input_file.close();
    return;
};
//End: read_matrix

// Gets the degeneration for a specified orbital.
// ener      : a vector containing MO energies.             (in)
// n_orb_max : the total number of MOs.                     (in)
// n_orb     : the desired MO to check for degeneracy.      (in)
// n_deg     : the degeneration number.                     (out)
// deg_MO    : a vector with the index for degenerated MOs. (out)
void get_degeneration(double *ener, int n_orb_max, int n_orb,
                      int *n_deg  , int *deg_MO)
{
    int    count, count_up, count_down;
    double criterium, ratio;

    // If there's only one MO in the system, performs no calculation.
    if (n_orb_max < 2) 
    {
        *n_deg    = 1;
        deg_MO[0] = 1;
    };

    * n_deg = 0; ratio = 0;
    criterium = 0.000005;

    if (n_orb > 1)
    {
        count_up   = n_orb  ;
        count_down = n_orb-1;
    } else {
        count_up   = 2;
        count_down = 1;
    };

    // Starts comparing bottom-up.
    for (count = count_up; count < n_orb_max+1; count += 1)
    {
        ratio = fabs(2 * (ener[n_orb] - ener[count]) 
                       / (ener[n_orb] + ener[count]));
        if (ratio < criterium)
        {
            *n_deg += 1;
            deg_MO[*n_deg - 1] = count;
        } else {
            break;
        };
    };

    // Starts comparing top-down.
    for (count = count_down; count > 0; count -= 1)
    {
        ratio = fabs(2 * (ener[n_orb] - ener[count])
                       / (ener[n_orb] + ener[count]));
        if (ratio < criterium)
        {
            *n_deg += 1;
            deg_MO[*n_deg - 1] = count;
        } else {
            break;
        };
    };
    return;
}; // get_degeneration

// Calculates the shape factor for a specified orbital.
// n_functs    : the number of total MOs.       (in)
// n_MO        : the MO to calculate the shape. (in)
// at_of_orb   : the basis function's atom.     (in)
// s_matrix    : the overlap matrix.            (in)
// coeff_mat   : the MO coefficient matrix.     (in)
// degen       : the degeneration of the MO.    (in)
// degen_MO    : the other generated MOs.       (in)
// shape_factor: the MO's shape factor.         (out)
void get_shape_factor(int n_functs, int * at_of_orb,
                      double * s_matrix, double * coeff_mat, 
                      int degen, int * degen_MO, double * shape_factor)
{
    int index_x, index_y, i_deg,   // Counters.
        index_a, index_b, index_s; // Matrix indexes.
    double dummy;

    for (index_x = 0; index_x < n_functs; index_x +=1){
        shape_factor[at_of_orb[index_x]-1] = 0;
    };    

    for (index_x = 0; index_x < n_functs; index_x += 1){
        for (index_y = 0; index_y < n_functs; index_y += 1){
            for (i_deg = 0; i_deg < degen; i_deg +=1){
                
                index_a = index_x + n_functs*(degen_MO[i_deg]);
                index_b = index_y + n_functs*(degen_MO[i_deg]);
                index_s = index_y + n_functs*(index_x);

                dummy = coeff_mat[index_a] * coeff_mat[index_b]
                                           * s_matrix[index_s];
                shape_factor[at_of_orb[index_x]-1] += dummy / degen;
            }
        };

    };
    return;

}; // get_shape_factor


// Calculates the shape factor for a specified orbital.
// shape_AH, shape_AL : Alpha shape factors.    (in)
// shape_BH, shape_BL : Beta shape factors.     (in)
// n_atoms            : The number of atoms.    (in)
// fukuiXX            : The fukui functions.    (out)
void get_fukui(double * shape_AH, double * shape_AL, double * shape_BH, 
               double * shape_BL, double * fukui_nn, double * fukui_ss, 
               double * fukui_ns, double * fukui_sn, int n_atoms)
{
    int    index;
    
    for (index = 0; index < n_atoms; index += 1)
    {
        fukui_nn[index + n_atoms]   = (shape_AH[index] + shape_BH[index]) / 2;
        fukui_nn[index + 2*n_atoms] = (shape_AL[index] + shape_BL[index]) / 2;
        fukui_nn[index]             = (fukui_nn[index + 2*n_atoms] +
                                       fukui_nn[index + n_atoms]) / 2;

        fukui_ss[index + n_atoms]   = (shape_AH[index] + shape_BL[index]) / 2;
        fukui_ss[index + 2*n_atoms] = (shape_AL[index] + shape_BH[index]) / 2;
        fukui_ss[index]             = (fukui_ss[index + 2*n_atoms] +
                                       fukui_ss[index + n_atoms]) / 2;

        fukui_ns[index + n_atoms]   = (shape_AH[index] - shape_BL[index]) / 2;
        fukui_ns[index + 2*n_atoms] = (shape_AL[index] - shape_BH[index]) / 2;
        fukui_ns[index]             = (fukui_ns[index + 2*n_atoms] +
                                       fukui_ns[index + n_atoms]) / 2;

        fukui_sn[index + n_atoms]   = (shape_AH[index] - shape_BH[index]) / 2;
        fukui_sn[index + 2*n_atoms] = (shape_AL[index] - shape_BL[index]) / 2;
        fukui_sn[index]             = (fukui_sn[index + 2*n_atoms] +
                                       fukui_sn[index + n_atoms]) / 2;
    };

}; // get_fukui

// Gets the system's global softness.
// enAH, enAL: the alpha HOMO and LUMO energies. (in)
// enBH, enBL: the beta HOMO and LUMO energies.  (in)
double get_softness(double enAH, double enAL, double enBH, double enBL)
{
    double softness; // The global softness.

    softness = 4 / (enAH + enBH - enAL - enBL);

    return softness;
};

// Formats a number into string, and adds a blank if positive.
std::string format_d(double number)
{
   std::string numb_str;
   
   if (number > 0){
        numb_str = " " + std::to_string(number);
   } else {
        numb_str = std::to_string(number);
   };

   return numb_str;
};

// Formats a number into string, and adds a blank if lesser than 10.
std::string format_i(int number)
{
   std::string numb_str;
   
   if (number < 10){
        numb_str = " " + std::to_string(number);
   } else {
        numb_str = std::to_string(number);
   };

   return numb_str;
};

// Main program.
int main(int argc, char *argv[])
{
    int    index_x, index_y, index_mat,  // Indexes for matrixes.
           n_orbs, n_elec_a, n_elec_b,   // Number of MOs and alpha/beta e-.
           n_atoms,                      // Number of atoms in molecule.
           n_deg_AH, n_deg_AL,           // Degeneration for alpha HOMO/LUMO.
           n_deg_BH, n_deg_BL;           // Degeneration for beta HOMO/LUMO.
    double softness;                     // The global softness.
    int    * atom_of_orb,                // The atom for each basis function.
           * deg_AH_MO, * deg_AL_MO ,    // Degenerated alpha HO/LU MOs.
           * deg_BH_MO, * deg_BL_MO ;    // Degenerated beta HO/LU MOs.
    double * energ_a, * energ_b,         // Alpha and beta energies.
           * coeff_a, * coeff_b,         // Alpha and beta MO coefficients.
           * s_matrix,                   // Overlap matrix.
           * shape_AH, * shape_AL,       // Alpha HO/LU shape factors.
           * shape_BH, * shape_BL,       // Beta HO/LU shape factors.
           * fukui_nn, * fukui_ss,       // The fukui NN and SS functions.
           * fukui_ns, * fukui_sn;       // The fukui NS and SN functions.

    std::cout << "\n########################################\n";
    std::cout << "#      FUKUI FUNCTION CALCULATION      #\n";
    std::cout << "########################################\n\n";

    n_orbs = 0; n_elec_a = 0; n_elec_b = 0;

    if (argc != 2)
    { 
        std::cout << "ERROR: Input filename missing from arguments.\n";
    } else {
        // Reads the number of orbital basis functions and the number
        // of alpha and beta electrons.
        read_n_orbs(argv[1], &n_orbs, &n_elec_a, &n_elec_b);
        std::cout << "Number of basis functions: " << n_orbs << ".\n";
        std::cout << "Number of Alpha electrons: " << n_elec_a << ".\n";
        std::cout << "Number of Beta  electrons: " << n_elec_b << ".\n";
      
        // Reads the energies, and overlap/MO-coefficient matrixes.
        coeff_a  = new double[n_orbs*n_orbs]; energ_a  = new double[n_orbs];
        coeff_b  = new double[n_orbs*n_orbs]; energ_b  = new double[n_orbs];
        s_matrix = new double[n_orbs*n_orbs]; atom_of_orb = new int[n_orbs];
        read_matrix(argv[1], n_orbs, s_matrix, atom_of_orb,
                    energ_a, energ_b, coeff_a, coeff_b);

        // Gets the number of atoms.
        n_atoms = atom_of_orb[n_orbs-1];
        std::cout << "Number of atoms:           " << n_atoms << ".\n\n";

        // Gets the degeneration of HOMO and LUMO orbitals.
        deg_AH_MO = new int[n_orbs]; deg_AL_MO = new int[n_orbs];
        deg_BH_MO = new int[n_orbs]; deg_BL_MO = new int[n_orbs];
        get_degeneration(energ_a, n_orbs, n_elec_a - 1, &n_deg_AH, deg_AH_MO);
        get_degeneration(energ_a, n_orbs, n_elec_a    , &n_deg_AL, deg_AL_MO);
        get_degeneration(energ_b, n_orbs, n_elec_b - 1, &n_deg_BH, deg_BH_MO);
        get_degeneration(energ_b, n_orbs, n_elec_b    , &n_deg_BL, deg_BL_MO);

        std::cout << "HOMO Alpha: " << format_d(energ_a[n_elec_a - 1])
                  << " Hartree. Degeneration: " << n_deg_AH <<  ". \n";
        std::cout << "LUMO Alpha: " << format_d(energ_a[n_elec_a])
                 << " Hartree. Degeneration: " << n_deg_AL <<  ". \n";
        std::cout << "HOMO Beta:  " << format_d(energ_b[n_elec_b - 1])
                 << " Hartree. Degeneration: " << n_deg_BH <<  ". \n";
        std::cout << "LUMO Beta:  " << format_d(energ_b[n_elec_b])
                 << " Hartree. Degeneration: " << n_deg_BL <<  ". \n";

        // Calculates the HOMO and LUMO orbital shape factors.
        shape_AH = new double[n_atoms]; shape_AL = new double[n_atoms];
        shape_BH = new double[n_atoms]; shape_BL = new double[n_atoms];

        get_shape_factor(n_orbs, atom_of_orb, s_matrix, coeff_a, n_deg_AH,
                         deg_AH_MO, shape_AH);
        get_shape_factor(n_orbs, atom_of_orb, s_matrix, coeff_a, n_deg_AL,
                         deg_AL_MO, shape_AL);
        get_shape_factor(n_orbs, atom_of_orb, s_matrix, coeff_b, n_deg_BH,
                         deg_BH_MO, shape_BH);
        get_shape_factor(n_orbs, atom_of_orb, s_matrix, coeff_b, n_deg_BL,
                         deg_BL_MO, shape_BL);

        // Calculates the Fukui functions.
        fukui_nn = new double[3*n_atoms]; fukui_ns = new double[3*n_atoms]; 
        fukui_ss = new double[3*n_atoms]; fukui_sn = new double[3*n_atoms]; 
        get_fukui(shape_AH, shape_AL, shape_BH, shape_BL, 
                  fukui_nn, fukui_ss, fukui_ns, fukui_sn, n_atoms);

        // Gets the molecule's global softness.
        softness = get_softness(energ_a[n_elec_a-1], energ_a[n_elec_a], 
                                energ_b[n_elec_b-1], energ_b[n_elec_b]);
        std::cout << "\nGlobal softness: " << softness << ".\n";
        
        std::cout << "\nAtom | FukuiNN- | FukuiNN+ | FukuiNN0 | Local softness \n";
        for (index_x = 0; index_x < n_atoms; index_x +=1)
        {
            std::cout << " " << format_i(index_x)      << "  | " 
                      << fukui_nn[index_x + n_atoms]   << "  " 
                      << fukui_nn[index_x + 2*n_atoms] << "  "
                      << fukui_nn[index_x]             << "  "
                      << fukui_nn[index_x]*softness    << "\n";
        };

        std::cout << "\nAtom | FukuiSS- | FukuiSS+ | FukuiSS0 | Local softness \n";
        for (index_x = 0; index_x < n_atoms; index_x +=1)
        {
            std::cout << " " << format_i(index_x)      << "  | " 
                      << fukui_ss[index_x + n_atoms]   << "  " 
                      << fukui_ss[index_x + 2*n_atoms] << "  "
                      << fukui_ss[index_x]             << "  "
                      << fukui_ss[index_x]*softness    << "\n";
        };

        std::cout << "\nAtom | FukuiNS- | FukuiNS+ | FukuiNS0 \n";
        for (index_x = 0; index_x < n_atoms; index_x +=1)
        {
            std::cout << " " << format_i(index_x)      << "  | " 
                      << fukui_ns[index_x + n_atoms]   << "  " 
                      << fukui_ns[index_x + 2*n_atoms] << "  "
                      << fukui_ns[index_x]             << "\n";
        };

        std::cout << "\nAtom | FukuiSN- | FukuiSN+ | FukuiSN0 \n";
        for (index_x = 0; index_x < n_atoms; index_x +=1)
        {
            std::cout << " " << format_i(index_x)      << "  | " 
                      << fukui_sn[index_x + n_atoms]   << "  " 
                      << fukui_sn[index_x + 2*n_atoms] << "  "
                      << fukui_sn[index_x]             << "\n";
        };
        
    };
    std::cout << "\n";
};


