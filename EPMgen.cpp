# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <algorithm>
#include <map>
#include <set>
#include <tuple>
#include <utility>   // for std::pair
# include <random>
# include <ctime>
using namespace std;

struct atom
{
    long int id;                        // atom id, globally unique
    atom(long int _id) : id(_id) {}     // Constructor to initialize atom id with _id
    int type;                           // atom type  1 for CH, 2 for CH2, 3 for CH3
    int monomerid;                      // monomer index in chain
    int polymerid;                      // polymer index
    vector<long int> neighbors;         // neighboring atom unique id
    int x, y, z;                        // wrapped coordinates
    int rx,ry,rz;                       // unwrapped coordinates, real positions

    bool operator<(const atom& other) const { return id < other.id; } 
                                        // Operator to compare atom ids for sorting for set and map

    bool operator==(const atom& other) const { return id == other.id; } 
                                        // Operator to check equality of atom ids
};
struct monomer
{
    vector <atom> atoms;                // atoms in the monomer
    std::map<long int, std::vector<long int>> connections; 
                                        // connections to other monomers, 
                                        // key is atom id, value is vector of connected atom ids

    int monomerType;                    // type of monomers in the polymer 1:monomer_C2, 2:monomer_C3
    monomer(int _monomerType) : monomerType(_monomerType) {} 
                                        // Constructor to initialize monomer type

    int atomNum = 0;                    // number of atoms in the monomer
    bool head = 0;                      // 1 for head, 0 for non-head monomer
    bool tail = 0;                      // 1 for tail, 0 for non-tail monomer

    void addConnection(long int atom1_id, long int atom2_id)
    {
        // Add a connection between two atoms in the monomer
        auto& list1 = connections[atom1_id];
        if (std::find(list1.begin(), list1.end(), atom2_id) == list1.end()) {
            list1.push_back(atom2_id);
        }
        auto& list2 = connections[atom2_id];
        if (std::find(list2.begin(), list2.end(), atom1_id) == list2.end()) {
            list2.push_back(atom1_id);
        }
    }
    std::vector<long int> getNeighbors(long int atom_id) const
    {
        // Get the neighbors of a given atom in the monomer
        if (connections.count(atom_id)) {
            return connections.at(atom_id);
        }
        return {}; // Return an empty vector if the atom_id is not found in connections
    }
};
struct polymer
{
    vector <monomer> monomers;              // monomers in the polymer
    std::map<long int, std::vector<long int>> monomerConnections; 
                                            // monomer connections between atoms in different monomers
                                            
    long int atomNum = 0;                        // number of atoms in the polymer
    int monomerNum = 0;                     // number of monomers in the polymer
    long int bondNum = 0;                   // number of bonds in the polymer
    long int angleNum = 0;                  // number of angles in the polymer
    long int dihedralNum = 0;               // number of dihedrals in the polymer
    long int improperNum = 0;               // number of impropers in the polymer
    void addMonomerConnection(long int atom1_id, long int atom2_id)
    {
        // Add a monomer connection between two atoms in different monomers
        auto& list1 = monomerConnections[atom1_id];
        if (std::find(list1.begin(), list1.end(), atom2_id) == list1.end()) {
            list1.push_back(atom2_id);
        }
        auto& list2 = monomerConnections[atom2_id];
        if (std::find(list2.begin(), list2.end(), atom1_id) == list2.end()) {
            list2.push_back(atom1_id);
        }
    }
    std::vector<long int> getMonomerNeighbors(long int atom_id) const
    {
        // Get the neighbors of a given atom in the polymer
        if (monomerConnections.count(atom_id)) {
            return monomerConnections.at(atom_id);
        }
        return {}; // Return an empty vector if the atom_id is not found in connections
    }
    std::set<long int> getAllatomIDs() const
    {
        // Helper to get all atom IDs in the polymer
        std::set<long int> atom_ids;
        for (const auto& monomer : monomers) {
            for (const auto& atom : monomer.atoms) {
                atom_ids.insert(atom.id);
            }
        }
        return atom_ids;
    }
};
struct list2
{
    long int id;
    long int type;
    long int atomid1;
    long int atomid2;
};
struct list3
{
    long int id;
    long int type;
    long int atomid1;
    long int atomid2;
    long int atomid3;
};
struct list4_line
{
    long int id;
    long int type;
    long int atomid1;
    long int atomid2;
    long int atomid3;
    long int atomid4;
};
struct list4_star
{
    long int id;
    long int type;
    long int centerid;
    long int atomid1;
    long int atomid2;
    long int atomid3;
};
vector<list2> bondList;          // Define a vector of list2 globally
vector<list3> angleList;         // Define a vector of list3 globally
vector<list4_line> dihedralList; // Define a vector of list4_line globally
vector<list4_star> improperList; // Define a vector of list4_star globally
struct boxLattice
{
    bool fcc;                   // 1: free, 0: occupied
    bool free_mol = 1;
    bool free_all = 1;
};
vector<polymer> epm;                   // Define a vector of polymer globally (will resize later)
vector<vector<vector<boxLattice>>> box;     // Define a 3D vector of boxLattice globally (will resize later)
void molFree(int n)
{
    // init the mol occupancy, free the molecule/polymer chain occupancy in box
    for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) for(int k = 0; k < n; k++) box[i][j][k].free_mol = 1;
}
ofstream outfile("monomer_types.txt");
void find2connections()
{
    // Find all 2-connection atoms in the polymers
    outfile << "--- Connected Two Atoms ---" << std::endl;
    long int bondid = 1; // Initialize bond id
    for (int i = 0; i < epm.size(); ++i) {
        const polymer& p = epm[i];
        outfile << "Polymer " << i << ":" << std::endl;
        std::set<std::pair<long int, long int>> reported_pairs; // To avoid duplicate reporting

        for (const auto& entry : p.monomerConnections) {
            long int atom1_id = entry.first;
            for (long int atom2_id : entry.second) { // A - B
                if (atom1_id < atom2_id) { // Ensure unique pairs (A, B) and (B, A) are not reported twice
                    if (reported_pairs.find({atom1_id, atom2_id}) == reported_pairs.end()) {
                        outfile << "  " << atom1_id << " - " << atom2_id << std::endl;
                        epm[i].bondNum++; // Increment bond count for the polymer
                        // Add to the bond list
                        bondList.push_back({bondid++, 1, atom1_id, atom2_id}); // type 1 for all bonds
                        reported_pairs.insert({atom1_id, atom2_id});
                    }
                }
            }
        }
    }
    long int totalBondNum = 0;
    for (const auto& polymer : epm) {
        totalBondNum += polymer.bondNum; // Sum up the bond counts from all polymers
    }
    outfile << "Total number of bonds: " << totalBondNum << std::endl; // Assuming all polymers have the same bond count
}
void find3connections()
{
    // find all 3-connection atoms in the polymers
    outfile << "--- Connected Three Atoms ---" << std::endl;
    long int angleid = 1; // Initialize angle id
    for (int i = 0; i < epm.size(); ++i) {
        const polymer& p = epm[i];
        outfile << "Polymer " << i << ":" << std::endl;
        std::set<std::vector<long int>> reported_triplets; // To avoid duplicate reporting
        for (const auto& entry : p.monomerConnections) {
            long int atom0_id = entry.first;
            for (long int atom1_id : entry.second) { // A - B
                for (long int atom2_id : p.getMonomerNeighbors(atom1_id)) { // B - C
                    if (atom0_id == atom2_id || atom1_id == atom2_id) continue; // Skip self-connections and duplicates
                    std::vector<long int> sorted_triplet = {atom0_id, atom1_id, atom2_id};
                    std::sort(sorted_triplet.begin(), sorted_triplet.end()); // Sort to ensure unique representation
                    if (reported_triplets.find(sorted_triplet) == reported_triplets.end()) {
                        outfile << "  " << atom0_id << " - " << atom1_id << " - " << atom2_id << std::endl;
                        epm[i].angleNum++; // Increment angle count for the polymer
                        // Add to the angle list
                        int atom1_type = -1;
                        // Search all monomers in the polymer for atom1_id to get the atom type
                        for (const auto& mon : epm[i].monomers) {
                            for (const auto& at : mon.atoms) {
                                if (at.id == atom1_id) {
                                    atom1_type = at.type;
                                    break;
                                }
                            }
                            if (atom1_type != -1) break;
                        }
                        if (atom1_type == 1)        // CH
                        {
                            angleList.push_back({angleid++, 1, atom0_id, atom1_id, atom2_id}); // type 1 for CH angles
                        } 
                        else if (atom1_type == 2)   // CH2
                        {
                            angleList.push_back({angleid++, 2, atom0_id, atom1_id, atom2_id}); // type 2 for CH2 angles
                        }
                        reported_triplets.insert(sorted_triplet);
                    }
                }
            }
        }
    }
    long int totalAngleNum = 0;
    for (const auto& polymer : epm) {
        totalAngleNum += polymer.angleNum; // Sum up the angle counts from all polymers
    }
    outfile << "Total number of angles: " << totalAngleNum << std::endl; // Assuming all polymers have the same angle count
}
void find4connections()
{
    // find all 4-connection atoms in the polymers
    outfile << "--- Connected Four Atoms in line ---" << std::endl;
    long int dihedralid = 1; // Initialize dihedral id
    for (int i = 0; i < epm.size(); ++i) {
        const polymer& p = epm[i];
        outfile << "Polymer " << i << ":" << std::endl;
        std::set<std::vector<long int>> reported_quadruplets; // To avoid duplicate reporting
        for (const auto& entry : p.monomerConnections) {
            long int atom0_id = entry.first;
            for (long int atom1_id : entry.second) { // A - B
                for (long int atom2_id : p.getMonomerNeighbors(atom1_id)) { // B - C
                    if (atom0_id == atom2_id) continue; // Skip self-connections
                    for (long int atom3_id : p.getMonomerNeighbors(atom2_id)) { // C - D
                        if (atom3_id == atom0_id || atom3_id == atom1_id) continue; // Skip duplicates
                        std::vector<long int> sorted_quadruplet = {atom0_id, atom1_id, atom2_id, atom3_id};
                        std::sort(sorted_quadruplet.begin(), sorted_quadruplet.end()); // Sort to ensure unique representation
                        if (reported_quadruplets.find(sorted_quadruplet) == reported_quadruplets.end()) {
                            outfile << "  " << atom0_id << " - " << atom1_id << " - " << atom2_id << " - " << atom3_id << std::endl;
                            epm[i].dihedralNum++; // Increment dihedral count for the polymer
                            // Add to the dihedral list
                            dihedralList.push_back({dihedralid++, 1, atom0_id, atom1_id, atom2_id, atom3_id}); // type 1 for all dihedrals
                            reported_quadruplets.insert(sorted_quadruplet);
                        }
                    }
                }
            }
        }
    }
    long int totalDihedralNum = 0;
    for (const auto& polymer : epm) {
        totalDihedralNum += polymer.dihedralNum; // Sum up the dihedral counts from all polymers
    }
    outfile << "Total number of dihedrals: " << totalDihedralNum << std::endl; // Assuming all polymers have the same dihedral count

    outfile << "--- Connected Four Atoms in star ---" << std::endl;
    long int improperid = 1; // Initialize improper id
    for (int i = 0; i < epm.size(); ++i) {
        const polymer& p = epm[i];
        outfile << "Polymer " << i << ":" << std::endl;
        std::set<std::vector<long int>> reported_stars; // To avoid duplicate reporting
        for (const auto& entry : p.monomerConnections) {
            long int center_id = entry.first;
            const auto& neighbors = entry.second;
            if (neighbors.size() < 3) continue; // Need at least 3 neighbors for a star
            // Generate all unique triplets of neighbors
            for (size_t a = 0; a < neighbors.size(); ++a) {
                for (size_t b = a + 1; b < neighbors.size(); ++b) {
                    for (size_t c = b + 1; c < neighbors.size(); ++c) {
                        std::vector<long int> star = {center_id, neighbors[a], neighbors[b], neighbors[c]};
                        std::sort(star.begin(), star.end());
                        if (reported_stars.find(star) == reported_stars.end()) {
                            outfile << "  " << center_id << " - (" << neighbors[a] << ", " << neighbors[b] << ", " << neighbors[c] << ")" << std::endl;
                            epm[i].improperNum++; // Increment improper count for the polymer
                            // Add to the improper list
                            improperList.push_back({improperid++, 1, center_id, neighbors[a], neighbors[b], neighbors[c]}); // type 1 for all impropers
                            reported_stars.insert(star);
                        }
                    }
                }
            }
        }
    }
    long int totalImproperNum = 0;
    for (const auto& polymer : epm) {
        totalImproperNum += polymer.improperNum; // Sum up the improper counts from all polymers
    }
    outfile << "Total number of impropers: " << totalImproperNum << std::endl; // Assuming all polymers have the same improper count
}

int main()
{
    long int totAtomNum = 0;
    srand(static_cast<unsigned int>(time(nullptr))); // Seed the random number generator
    int chainNum = 100;
    int chainLen = 250;
    double density = 0.05;
    epm.resize(chainNum); // Resize the vector to hold 'chainNum' polymers
    double m2rate = 0.3;
    for (int i = 0; i < chainNum; i++)
    {
        epm[i].monomerNum = chainLen;
        for (int j = 0; j < chainLen; j++)
        {
            double k = static_cast<double>(rand()) / RAND_MAX;
            if (k < m2rate)
            {
                epm[i].monomers.push_back(monomer(1)); // Create a monomer of type 1 (monomer_C2)
                epm[i].atomNum += 2;
                totAtomNum += 2;
            } 
            else
            {
                epm[i].monomers.push_back(monomer(2)); // Create a monomer of type 2 (monomer_C3)
                // which means epm[i].monomers.monomerType[j] = 2;
                epm[i].atomNum += 3;
                totAtomNum += 3;
            } 

        }
    }

    outfile << "Total number of atoms: " << totAtomNum << std::endl;
    long int atomId = 1;
    for (int i = 0; i < chainNum; i++)
    {
        for (int j = 0; j < epm[i].monomerNum; j++)
        {
            if (epm[i].monomers[j].monomerType == 1)                     // monomer_C2: -CH2-CH2-
            {
                epm[i].monomers[j].atoms.push_back(atom(atomId++));
                epm[i].monomers[j].atoms.push_back(atom(atomId++));
                if (j == 0) 
                {
                    epm[i].monomers[j].head = 1;                // first monomer is head
                    epm[i].monomers[j].atoms[0].type = 3;       // CH3
                    epm[i].monomers[j].atoms[1].type = 2;       // CH2
                }
                else if (j == epm[i].monomerNum - 1) 
                {
                    epm[i].monomers[j].tail = 1;                // last monomer is tail
                    epm[i].monomers[j].atoms[0].type = 2;       // CH2
                    epm[i].monomers[j].atoms[1].type = 3;       // CH3
                }
                else
                {
                    epm[i].monomers[j].head = 0;                // not head
                    epm[i].monomers[j].atoms[0].type = 2;       // CH2
                    epm[i].monomers[j].atoms[1].type = 2;       // CH2
                }
                epm[i].monomers[j].atoms[0].monomerid = j;      // Set monomer id 
                epm[i].monomers[j].atoms[0].polymerid = i;      // Set polymer id
                epm[i].monomers[j].atoms[1].monomerid = j;      // Set monomer id
                epm[i].monomers[j].atoms[1].polymerid = i;      // Set polymer id
                epm[i].monomers[j].atomNum = 2;                 // Set atom number in the monomer

                // Add connections between atoms in the monomer
                epm[i].monomers[j].addConnection(epm[i].monomers[j].atoms[0].id, epm[i].monomers[j].atoms[1].id);
            } 
            else if (epm[i].monomers[j].monomerType  == 2)               // monomer_C3: -CH2-CH(CH3)-
            {
                epm[i].monomers[j].atoms.push_back(atom(atomId++));
                epm[i].monomers[j].atoms.push_back(atom(atomId++));
                epm[i].monomers[j].atoms.push_back(atom(atomId++));
                if (j == 0) 
                {
                    epm[i].monomers[j].head = 1;                // first monomer is head
                    epm[i].monomers[j].atoms[0].type = 3;       // CH3
                    epm[i].monomers[j].atoms[1].type = 1;       // CH
                    epm[i].monomers[j].atoms[2].type = 3;       // CH3
                }
                else if (j == epm[i].monomerNum - 1) 
                {
                    epm[i].monomers[j].tail = 1;                // last monomer is tail
                    epm[i].monomers[j].atoms[0].type = 2;       // CH2
                    epm[i].monomers[j].atoms[1].type = 2;       // CH2
                    epm[i].monomers[j].atoms[2].type = 3;       // CH3
                }
                else
                {
                    epm[i].monomers[j].atoms[0].type = 2;       // CH2
                    epm[i].monomers[j].atoms[1].type = 1;       // CH
                    epm[i].monomers[j].atoms[2].type = 3;       // CH3
                }
                for (int k = 0; k < 3; k++)
                {
                    epm[i].monomers[j].atoms[k].monomerid = j;      // Set monomer id
                    epm[i].monomers[j].atoms[k].polymerid = i;      // Set polymer id
                }
                epm[i].monomers[j].atomNum = 3;                 // Set atom number in the monomer

                // Add connections between atoms in the monomer
                epm[i].monomers[j].addConnection(epm[i].monomers[j].atoms[0].id, epm[i].monomers[j].atoms[1].id);
                epm[i].monomers[j].addConnection(epm[i].monomers[j].atoms[1].id, epm[i].monomers[j].atoms[2].id);
            }
        }
        
        // add connections between monomers in a chain
        for (int j = 0; j < epm[i].monomerNum - 1; j++)
        {
            epm[i].addMonomerConnection(epm[i].monomers[j].atoms[1].id, epm[i].monomers[j + 1].atoms[0].id);
        }
        for (const auto& monomer : epm[i].monomers) {
            for (const auto& entry : monomer.connections) {
                long int atom0_id = entry.first;
                for (int atom1_id : entry.second) {
                    // Add connections between atoms in different monomers
                    epm[i].addMonomerConnection(atom0_id, atom1_id);
                }
            }
        }

    }



    // let's go self-avoid random walk
    long int boxLen = 2;
    long long boxVol = boxLen * boxLen * boxLen * density;
    while (boxVol <= totAtomNum)
    {
        boxLen += 2;
        boxVol = boxLen * boxLen * boxLen * density;
    } 
    box.resize(boxLen);                     // Resize the 3D vector 'box' to boxLen in each direction
    for (int i = 0; i < boxLen; ++i)
    {
        box[i].resize(boxLen);
        for (int j = 0; j < boxLen; ++j)
        {
            box[i][j].resize(boxLen);
        }
    }
    
    for (int i = 0; i < boxLen; i++)        // init the fcc box
    {
        for (int j = 0; j < boxLen; j++)
        {
            for (int k = 0; k < boxLen; k++)
            {
                if (k % 2 == 0)
                {
                    box[i][j][k].fcc = (i + j) % 2;
                    box[i][j][k].free_all = 1; // Initialize all boxes as free
                    box[i][j][k].free_mol = 1; // Initialize all boxes as free for molecules
                }
                else
                {
                    box[i][j][k].fcc = (i + j + 1) % 2;
                    box[i][j][k].free_all = 1; // Initialize all boxes as free
                    box[i][j][k].free_mol = 1; // Initialize all boxes as free for molecules
                }
            }
        }
    }

    // Assign atoms to boxes with self-avoid random walk

    // outfile << "Box length: " << boxLen << endl;
    // for (int i = 0; i < chainNum; i++)
    // {
    //     outfile << epm[i].atomNum << "      ";
    //     for (int j = 0; j < chainLen; j++) {
    //         outfile << epm[i].monomers[j].monomerType << "    ";
    //     }
    //     outfile << endl;
    // }

    for (int i = 0; i < chainNum; i++)                              // polymer chain loop
    {
        // outfile << "Placing polymer chain " << i + 1 << " with " << epm[i].monomerNum<< " monomers." << endl;
        molFree(boxLen);
        bool Trapped = 0;                                           // Flag to check if the polymer chain is trapped
        for (int j = 0; j < epm[i].monomerNum; j++)                 // monomer loop
        {
            // outfile << "Placing monomer " << j + 1 << " of type " << epm[i].monomers[j].monomerType << endl;
            for (int k = 0; k < epm[i].monomers[j].atomNum; k++)    // atom loop in the monomer
            {
                if (j == 0 && k == 0)                       // random walk for the first atom
                {
                    int wx, wy, wz;
                    do {
                        wx = rand() % boxLen;
                        wy = rand() % boxLen;
                        wz = rand() % boxLen;
                    } while (!box[wx][wy][wz].fcc || !box[wx][wy][wz].free_all);
                    epm[i].monomers[j].atoms[k].x = wx;
                    epm[i].monomers[j].atoms[k].y = wy;
                    epm[i].monomers[j].atoms[k].z = wz;
                    epm[i].monomers[j].atoms[k].rx = wx;
                    epm[i].monomers[j].atoms[k].ry = wy;
                    epm[i].monomers[j].atoms[k].rz = wz;
                    box[wx][wy][wz].free_mol = 0;
                    // outfile << "Placing first atom at: (" << wx << ", " << wy << ", " << wz << ")" << endl;
                }
                else // random walk for the rest atoms
                {
                    // Find the previous atom's position
                    int prev_j = j, prev_k = k - 1;
                    if (prev_k < 0) { prev_j = j - 1; prev_k = 1; }
                    int px = epm[i].monomers[prev_j].atoms[prev_k].x;
                    int py = epm[i].monomers[prev_j].atoms[prev_k].y;
                    int pz = epm[i].monomers[prev_j].atoms[prev_k].z;

                    // 12-connected neighbors for FCC lattice
                    vector<tuple<int, int, int>> directions = {
                        // Face-centered cubic (FCC) nearest neighbors
                        {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0},
                        {1, 0, 1}, {1, 0, -1}, {-1, 0, 1}, {-1, 0, -1},
                        {0, 1, 1}, {0, 1, -1}, {0, -1, 1}, {0, -1, -1}
                    };
                    vector<tuple<int, int, int, int, int, int>> candidates;        // find the unoccupied position for candidates
                    for (auto [dx, dy, dz] : directions) {
                        int nx = (int(px) + dx + boxLen) % boxLen;
                        int ny = (int(py) + dy + boxLen) % boxLen;
                        int nz = (int(pz) + dz + boxLen) % boxLen;
                        if (box[nx][ny][nz].fcc && box[nx][ny][nz].free_mol && box[nx][ny][nz].free_all)
                        {
                            candidates.push_back({nx, ny, nz, dx, dy, dz}); // Store both position and direction
                        }
                    }
                    if (candidates.empty()) {
                        Trapped = 1;
                        break; // No available position, trapped
                    }
                    // Randomly pick one available neighbor
                    int idx = rand() % candidates.size();
                    int dx, dy, dz;
                    int xpos, ypos, zpos;
                    std::tie(xpos, ypos, zpos, dx, dy, dz) = candidates[idx];
                    epm[i].monomers[j].atoms[k].x = xpos;
                    epm[i].monomers[j].atoms[k].y = ypos;
                    epm[i].monomers[j].atoms[k].z = zpos;
                    // Update unwrapped coordinates by adding direction to previous atom's rx, ry, rz
                    epm[i].monomers[j].atoms[k].rx = epm[i].monomers[prev_j].atoms[prev_k].rx + dx;
                    epm[i].monomers[j].atoms[k].ry = epm[i].monomers[prev_j].atoms[prev_k].ry + dy;
                    epm[i].monomers[j].atoms[k].rz = epm[i].monomers[prev_j].atoms[prev_k].rz + dz;
                    box[xpos][ypos][zpos].free_mol = 0;
                }
            }
            if (Trapped) break; // Break the monomer loop if trapped
        }
        if (Trapped)
        {
            i--; // Decrement i to retry this polymer chain
            continue; // Skip to the next iteration of the polymer chain loop
        }
        else
        {
            // After placing all atoms in the polymer chain, mark the box as occupied
            for (int j = 0; j < epm[i].monomerNum; j++)
            {
                for (int k = 0; k < epm[i].monomers[j].atomNum; k++)
                {
                    int x = epm[i].monomers[j].atoms[k].x;
                    int y = epm[i].monomers[j].atoms[k].y;
                    int z = epm[i].monomers[j].atoms[k].z;
                    box[x][y][z].free_all = 0; // Mark the box as occupied
                }
            }
        }
    }
    // find2connections();
    // outfile << std::endl;
    // find3connections();
    // outfile << std::endl;
    // find4connections();

    // Generate bonds, angles, dihedrals, and impropers
    find2connections();
    outfile << std::endl;
    find3connections();
    outfile << std::endl;
    find4connections();
    outfile << std::endl;

    // Output LAMMPS data file format
    ofstream lmpfile("EPM.data");
    lmpfile << "LAMMPS data file via EPM generator by Leo-LYY's code\n\n";
    long int totalBonds = bondList.size();
    long int totalAngles = angleList.size();
    long int totalDihedrals = dihedralList.size();
    long int totalImpropers = improperList.size();
    long int totalAtoms = totAtomNum;

    // Count types
    int atomTypes = 3;      // 1:CH, 2:CH2, 3:CH3
    int bondTypes = 1;
    int angleTypes = 2;
    int dihedralTypes = 1;
    int improperTypes = 1;

    // Box bounds (use boxLen for bounds)
    int rxmin = 0, rymin = 0, rzmin = 0;
    double sFactor = 1.08894; // 1.54/√2
    double rxmax = boxLen * sFactor, rymax = boxLen * sFactor, rzmax = boxLen * sFactor;

    lmpfile << totalAtoms << " atoms\n";
    lmpfile << atomTypes << " atom types\n";
    lmpfile << totalBonds << " bonds\n";
    lmpfile << bondTypes << " bond types\n";
    lmpfile << totalAngles << " angles\n";
    lmpfile << angleTypes << " angle types\n";
    lmpfile << totalDihedrals << " dihedrals\n";
    lmpfile << dihedralTypes << " dihedral types\n";
    lmpfile << totalImpropers << " impropers\n";
    lmpfile << improperTypes << " improper types\n\n";

    lmpfile << rxmin << " " << rxmax << " xlo xhi\n";
    lmpfile << rymin << " " << rymax << " ylo yhi\n";
    lmpfile << rzmin << " " << rzmax << " zlo zhi\n\n";

    // Atoms section
    lmpfile << "Atoms\n\n";
    // atomid molid type rx ry rz
    for (const auto& p : epm) {
        int molid = &p - &epm[0] + 1;
        for (const auto& m : p.monomers) {
            for (const auto& a : m.atoms) {
                double sx = a.rx * sFactor; // Scale the coordinates
                double sy = a.ry * sFactor;
                double sz = a.rz * sFactor;
                // Write atom data: id, molid, type, scaled coordinates
                lmpfile << a.id << " " << molid << " " << a.type << " "
                        << sx << " " << sy << " " << sz << "\n";
            }
        }
    }
    lmpfile << "\n";

    // Masses section
    lmpfile << "Masses\n\n";
    lmpfile << "1 13.019\n"; // Mass for CH
    lmpfile << "2 14.027\n"; // Mass for CH2
    lmpfile << "3 15.035\n"; // Mass for CH3
    lmpfile << "\n";

    // Bonds section
    lmpfile << "Bonds\n\n";
    for (const auto& b : bondList) {
        lmpfile << b.id << " " << b.type << " " << b.atomid1 << " " << b.atomid2 << "\n";
    }
    lmpfile << "\n";

    // Angles section
    lmpfile << "Angles\n\n";
    for (const auto& a : angleList) {
        lmpfile << a.id << " " << a.type << " " << a.atomid1 << " " << a.atomid2 << " " << a.atomid3 << "\n";
    }
    lmpfile << "\n";

    // Dihedrals section
    lmpfile << "Dihedrals\n\n";
    for (const auto& d : dihedralList) {
        lmpfile << d.id << " " << d.type << " " << d.atomid1 << " " << d.atomid2 << " " << d.atomid3 << " " << d.atomid4 << "\n";
    }
    lmpfile << "\n";

    // Impropers section
    lmpfile << "Impropers\n\n";
    for (const auto& im : improperList) {
        lmpfile << im.id << " " << im.type << " " << im.centerid << " " << im.atomid1 << " " << im.atomid2 << " " << im.atomid3 << "\n";
    }
    lmpfile << "\n";

    lmpfile.close();
    outfile.close();

}