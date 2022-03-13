#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <complex>
#include <chrono>
#include "Eigen/Dense"

#define pi 3.14159265359

double splitVal (std::string str){
    // value is the impedance with magnitude
    // mag is the magnitude of the impedance

    double value;
    char mag;
     
     mag = str[str.length() - 1];

     if (mag <= 57){
       value = std::stod(str);
       return value;
     }

     if (mag != 'g'){
         str.pop_back();
         value = std::stod(str);
     }
     // Special case for Meg
     else {
         str.pop_back();
         str.pop_back();
         str.pop_back();
         value = std::stod(str);
     }
    
    // To multiply according to magnitude
    if (mag == 'p'){
        value *= pow(10,-12);
    }
    else if (mag == 'n'){
        value *= pow(10,-9);
    }
    else if (mag == 'u'){
        value *= pow(10,-6);
    }
    else if (mag == 'm'){
        value *= pow(10,-3);
    }
    else if (mag == 'k'){
        value *= pow(10,3);
    }
    else if (mag == 'g'){
        value *= pow(10,6);
    }
    else if (mag == 'G'){
        value *= pow(10,9);
    }

    return value;
}

//checks whether the character is A
bool check_A (char str){
    if (str == 'A'){
        return true;
    }
    return false;
}

//decodes the netlist into a vector
std::vector<std::string> decodeNet (std::string string){
    std::vector <std::string> netList;
    std::string s;

    for (int j = 0; j < string.length(); j++){
        s.clear();
        while (string[j] != ' ' && j < string.length() && !check_A(string[j])){
            s += string[j];
            j++;
        }
        if (check_A(string[j])){
            while (j < string.length()){
                s += string[j];
                j++;
            }
        }
        netList.push_back(s);
    }
    
    return netList;
}

//separates AC(1 0) to 1 and 0 in a vector
std::vector<double> convertAmpPhase (std::string str){
    std::vector<double> output;
    std::string s;
    int i = 0;

    str.pop_back();
    str.erase(str.begin());
    str.erase(str.begin());
    str.erase(str.begin());
      
    while (str[i] != ' '){
        s += str[i];
        i++;
    }
    output.push_back(std::stoi(s));
    s.clear();

    while (i < str.length()){
        s += str[i];
        i++;
    }
    output.push_back(std::stoi(s));
      
    return output;
}

//returns the node number
int node (std::string str){
    int n;
    if (str == "0"){
        n = 0;
        return n;
    }
    else {
        str.erase(0,3);
        n = std::stoi(str);
    }

    return n;
}

//changes a vector into SPICE format
void addNetListResistor (std::vector<std::string>& netList, std::string node1, std::string node2, double res){
    netList.push_back("R");
    netList.push_back(node1);
    netList.push_back(node2);
    netList.push_back(std::to_string(res));
}

//changes a vector into SPICE format
void addNetListVS (std::vector<std::string>& netList, int node1, int node2, double volt){
    netList.push_back("V");
    if (node1 == 0){
        netList.push_back("0");
    }
    else{
        netList.push_back("N00" + std::to_string(node1));
    }
    if (node2 == 0){
        netList.push_back("0");
    }
    else{
        netList.push_back("N00" + std::to_string(node2));
    }
    netList.push_back(std::to_string(volt));
}

//changes a vector into SPICE format
void addNetListCS (std::vector<std::string>& netList, std::string node1, std::string node2, std::string node3, std::string node4, double trans){
    netList.push_back("G");
    netList.push_back(node1);
    netList.push_back(node2);
    netList.push_back(node3);
    netList.push_back(node4);
    netList.push_back(std::to_string(trans));
}


class netList {
    public:
        int n0;
        int n1;
};

class voltageSource: public netList {
    public:
        double amplitude;
        double phase;
        double dcVal;
        double acVal;
        int ind;
        std::complex<double> complexVal;
        std::vector<voltageSource> v;
    
    voltageSource add (std::vector<std::string> netList){
        voltageSource volt;
        if (netList[3][0] == 'A'){
            volt.ind = 1;
            std::vector<double> vec = convertAmpPhase(netList[3]);
            volt.complexVal = std::polar (vec[0], (vec[1] * pi / 180));
        }
        else{
            volt.ind = 0;
            volt.dcVal = std::stoi(netList[3]);
        }
        volt.n0 = node(netList[1]);
        volt.n1 = node(netList[2]);
       
        return volt;
    }
    
    void addFreq (double freq){
        double co = amplitude * cos((2 * pi * freq) + phase);
        acVal = co;
    }

};

class currentSource: public netList{
    public:
        double amplitude;
        double phase;
        double dcVal;
        std::vector <currentSource> c;
    
    currentSource add (std::vector<std::string> netList){
        currentSource curr;
        if (netList[3][0] == 'A'){
            std::vector<double> vec = convertAmpPhase(netList[3]);
            curr.amplitude = vec[0];
            curr.phase = vec[1];
        }
        else{
            double d = std::stod(netList[3]);
            curr.dcVal = d;
        }
        curr.n0 = node(netList[1]);
        curr.n1 = node(netList[2]);
       
        return curr;
    }
    
};

class resistor: public netList{
    public:
        double res;
    
    resistor add (std::vector<std::string> netList){
        resistor r;
        r.n0 = node(netList[1]);
        r.n1 = node(netList[2]);
        r.res = splitVal(netList[3]);
       
        return r;
    }
    
    std::complex<double> imp(){
        std::complex<double> impedance(res);
        
        return impedance;
    }
    
};

class capacitor: public netList{
    public:
        double cap;
    
    capacitor add (std::vector<std::string> netList){
        capacitor c;
        c.n0 = node(netList[1]);
        c.n1 = node(netList[2]);
        c.cap = splitVal(netList[3]);
        
        return c;
    }
    
    std::complex<double> imp(double freq){
        double omega = 2 * pi * freq;
        std::complex<double> impedance (0, - 1/(omega * cap));
        
        return impedance;
    }
};

class inductor: public netList{
    public:
        double induc;
    
    inductor add (std::vector<std::string> netList){
        inductor i;
        i.n0 = node(netList[1]);
        i.n1 = node(netList[2]);
        i.induc = splitVal(netList[3]);
        
        return i;
    }
    
    std::complex<double> imp(double freq){
        double omega = 2 * pi * freq;
        std::complex<double> impedance (0, omega * induc);
        
        return impedance;
    }
    
};

class diode: public netList{
    public:
        std::string model;

    diode add (std::vector<std::string> netList){
        diode d;
        d.n0 = node(netList[1]);
        d.n1 = node(netList[2]);
        
        return d;
    }
    
    void ssem (diode d, std::vector<std::string> netList, std::vector<resistor>& res){
        std::vector<std::string> netListr;
        resistor r;
        double id;
        double vd = 0.7;
        double is = 2.52 * pow(10,-9);
        double vt = 0.025;
        double n = 1.752;
    
        id = is * (exp(vd / (n * vt)) - 1);

        double rd = vt / id;
        
        std::string node_an = netList[1];
        std::string node_ca = netList[2];

        // for diode resistor
        addNetListResistor(netListr, node_an, node_ca, rd);
        res.push_back(r.add(netListr));

        netListr.clear();
    }

    void newton (double& vd, std::vector<resistor>& res){
        resistor r;
        std::vector<std::string> netListr;
        double is = 2.52 * pow(10,-9);
        double vt = 0.025;
        double n = 1.752;
        double id;
        double idd;
        double rd;
            
        id = is * (exp(vd / (n * vt)) - 1);
        idd = is * ((1/(n*vt)) * (exp(vd / (n * vt))) - 1);
        vd -= (id/idd);
        
        id = is * (exp(vd / (n * vt)) - 1);
        rd = vt / id;
        
        addNetListResistor(netListr, ("N00" + std::to_string(n0)), ("N00" + std::to_string(n1)), rd);
        res.push_back(r.add(netListr));

        netListr.clear();
    }
    
};

class controlSource: public netList{
    public:
        int n2;
        int n3;
        double transc;
    
    controlSource add (std::vector<std::string> netList){
        controlSource cs;
        cs.n0 = node(netList[1]);
        cs.n1 = node(netList[2]);
        cs.n2 = node(netList[3]);
        cs.n3 = node(netList[4]);
        cs.transc = splitVal(netList[5]);
        
        return cs;
    }
    
};

class bjt: public netList{
    public:
        int n2;
        std::string model;
    
    bjt add (std::vector<std::string> netList){
        bjt b;
        b.n0 = node(netList[1]);
        b.n1 = node(netList[2]);
        b.n2 = node(netList[3]);
        b.model = netList[4];
       
        return b;
    }
  
    void newtonBJT (double& vbe, double& vce, std::vector<resistor>& r, std::vector<controlSource>& g){
        double is = pow(10,-14);
        double vt = 0.025;
        double va = 100;
        double ic = pow(1,-3);
        double ib;
        double beta = 200;
        std::vector<std::string> netListr;
        std::vector<std::string> netListg;
        resistor res;
        controlSource cs;
       
        ic = is * exp(vbe/vt) * (1 + vce/va);
        ib = ic / beta;
        
        double rbe;
        double gm = ic / vbe;
        double ro;
        
        rbe = beta / gm;
        ro = va / ic;
        
        std::string node_c = "N00" + std::to_string(n0);
        std::string node_b = "N00" + std::to_string(n1);
        std::string node_e = "N00" + std::to_string(n2);

        //for rbe
        addNetListResistor(netListr, node_b, node_e, rbe);
        r.push_back(res.add(netListr));

        netListr.clear();
        
        //for ro
        addNetListResistor(netListr, node_c, node_e, ro);
        r.push_back(res.add(netListr));
        
        netListr.clear();

        //for current source
        addNetListCS(netListg, node_c, node_e, node_b, node_e, gm);
        g.push_back(cs.add(netListg));

        netListg.clear();
    }
    
    void ssem (bjt b, std::vector<controlSource>& g, std::vector<resistor>& r, std::vector<std::string> netList){
        double ic = 0.8;
        double vt = 0.025;
        double va = 100;
        double beta = 200;
        std::vector<std::string> netListr;
        std::vector<std::string> netListg;
        resistor res;
        controlSource cs;

        double rbe;
        double ro;
        double gm = ic / vt;
    
        rbe = beta / gm;
        ro = va / ic;
        
        std::string node_c = netList[1];
        std::string node_b = netList[2];
        std::string node_e = netList[3];
        
        //for rbe
        addNetListResistor(netListr, node_b, node_e, rbe);
        r.push_back(res.add(netListr));
        
        netListr.clear();
        
        //for ro
        addNetListResistor(netListr, node_c, node_e, ro);
        r.push_back(res.add(netListr));
        
        netListr.clear();
        
        //for current source
        addNetListCS(netListg, node_c, node_e, node_b, node_e, gm);
        g.push_back(cs.add(netListg));
        
        netListg.clear();
    }
};

class mosfet: public netList{
    public:
    int n2;
    std::string model;

    mosfet add (std::vector<std::string> netList){
        mosfet m;
        m.n0 = node(netList[1]);
        m.n1 = node(netList[2]);
        m.n2 = node(netList[3]);
        m.model = netList[4];
        return m;
    }
    
    void ssem (std::vector<std::string> netList, std::vector<resistor>& r, std::vector<controlSource>& g, std::vector<double> voltage){
        std::vector<std::string> netListr;
        std::vector<std::string> netListg;
        resistor rout;
        controlSource cs;
        
        int node_d = node(netList[1]);
        int node_g = node(netList[2]);
        int node_s = node(netList[3]);
        

        double vt = 0.025;
        double va = 150;
        double k = 10;
        
        if (model == "NMOS"){
            double vgs = voltage[node_g] - voltage[node_s];
            double id = k * pow((vgs - vt),2);
            
            double gm = 2 * sqrt(id * k);

            double ro = va / id;
            
            //for ro
            netListr.push_back("R");
            netListr.push_back(netList[1]);
            netListr.push_back(netList[3]);
            netListr.push_back(std::to_string(ro));
            r.push_back(rout.add(netListr));
            
            netListr.clear();
            
            //for control source
            netListg.push_back("G");
            netListg.push_back(netList[1]);
            netListg.push_back(netList[3]);
            netListg.push_back(netList[2]);
            netListg.push_back(netList[3]);
            netListg.push_back(std::to_string(gm * vgs));
            g.push_back(cs.add(netListg));
            
            netListg.clear();
        }
        
        else if (model == "PMOS"){
        
            double vsg = voltage[node_s] - voltage[node_g];
            double id = k * pow((vsg - vt),2);
            
            double gm = 2 * sqrt(id * k);

            double ro = va / id;
            
            //for ro
            netListr.push_back("R");
            netListr.push_back(netList[1]);
            netListr.push_back(netList[3]);
            netListr.push_back(std::to_string(ro));
            r.push_back(rout.add(netListr));
            
            netListr.clear();
            
            //for control source
            netListg.push_back("G");
            netListg.push_back(netList[3]);
            netListg.push_back(netList[1]);
            netListg.push_back(netList[3]);
            netListg.push_back(netList[2]);
            netListg.push_back(std::to_string(gm * vsg));
            g.push_back(cs.add(netListg));
            
            netListg.clear();
        }
    }
};

//Store each nodes into an int vector
std::vector<int> store_noder (std::vector<resistor> node){
    std::vector<int> nodesize;
    for (int i = 0; i < node.size(); i++){
        nodesize.push_back(node[i].n0);
        nodesize.push_back(node[i].n1);
    }
  
    return nodesize;
}

//Store each nodes into an int vector
std::vector<int> store_nodercl (std::vector<resistor> node1, std::vector<capacitor> node2, std::vector<inductor> node3){
    std::vector<int> nodesize;
    for (int i = 0; i < node1.size(); i++){
        nodesize.push_back(node1[i].n0);
        nodesize.push_back(node1[i].n1);
    }
    for (int i = 0; i <node2.size(); i++){
        nodesize.push_back(node2[i].n0);
        nodesize.push_back(node2[i].n1);
    }
    for (int i = 0; i < node3.size(); i++){
        nodesize.push_back(node3[i].n0);
        nodesize.push_back(node3[i].n1);
    }
    
    return nodesize;
}

void store_node (const std::vector<std::string> netList, std::vector<int>& nodesize){
    if (netList[0][0] != 'Q' || netList[0][0] != 'G' || netList[0][0] != 'M'){
        nodesize.push_back(node(netList[1]));
        nodesize.push_back(node(netList[2]));
    }
    else if (netList[0][0] == 'Q' || netList[0][0] == 'G'){
        nodesize.push_back(node(netList[1]));
        nodesize.push_back(node(netList[2]));
        nodesize.push_back(node(netList[3]));
    }
    else{
        nodesize.push_back(node(netList[1]));
        nodesize.push_back(node(netList[2]));
        nodesize.push_back(node(netList[3]));
        nodesize.push_back(node(netList[4]));
    }
}

// Find the highest number of the node
int max_node_size (std::vector<int> nodesize){
    int max = nodesize[0];
    for (int i = 0; i < nodesize.size(); i++){
        if (nodesize[i] > max){
            max = nodesize[i];
        }
    }
    
    return max;
}

std::vector <double> conduc_vector (std::vector<resistor> R){
    std::vector <double> conduc_val;
    for (int i = 0; i < R.size(); i++){
        conduc_val.push_back((1.0/R[i].res));
    }
    return conduc_val;
}

//a vector store all the impedance of the R , C and I
std::vector<std::complex<double>> overall_imp (std::vector<resistor> R, std::vector<capacitor> C, std::vector<inductor> I, double freq){
    std::vector<std::complex<double>> tmp;
    for (int i = 0; i < R.size(); i++){
        tmp.push_back(1.0 / R[i].imp());
    }
    for (int i = 0; i < C.size(); i++){
        tmp.push_back(1.0 / C[i].imp(freq));
    }
    for (int i = 0; i < I.size(); i++){
        tmp.push_back(1.0 / I[i].imp(freq));
    }
    
    return tmp;
}

// a function that store the resistor values in a vector when the resistor connected to the node n
std::vector<std::complex<double>> rcl_in_vec (int n, std::vector<int> node_list, std::vector<std::complex<double>> overall_imp, std::vector<int> node){
    std::vector<std::complex<double>> Node;
    for (int i = 0; i < node_list.size(); i++){
        if (node_list[i] == n){
            Node.push_back(overall_imp[i/2]);
        }
        if (n >= max_node_size(node) && (max_node_size(node_list) < max_node_size(node))){
            Node.push_back(0);
        }
    }
    
    int x = (int)Node.size();
    
    if (Node.size() < max_node_size(node)){
        for (int i = 0; i < (max_node_size(node)- x); i++){
            Node.push_back(0);
        }
    }
    
    return Node;
}

// a function that create a vector <vector <double>> and store the Node in first index and corresponding resistor value in second index
void Node_vec (const std::vector<int> &node_list, std::vector<std::complex<double>> overall_imp, std::vector<std::vector<std::complex<double>>> &rnode, std::vector<int> &node) {
    for (int i = 1; i <= max_node_size(node); i++){
        rnode.push_back(rcl_in_vec(i, node_list, overall_imp, node));
    }
}

//Finds the input in the inputed vector
int search_vector (std::complex<double> n, const std::vector<std::complex<double>> vin){
    for (int i=0; i < vin.size(); i++){
        if (vin[i] == n){
            return i;
        }
    }
    
  return -1;
}

//Finds whether there is an intersection between the 2 input vectors
void intersection_vector (const std::vector<std::complex<double>> &vin1, const std::vector<std::complex<double>> &vin2, std::vector<std::complex<double>> &vout){
    for (int i = 0; i < vin1.size(); i++){
        if (search_vector(vin1[i], vin2) != -1){
            vout.push_back(vin2[search_vector(vin1[i],vin2)]);
        }
    }
}

// function to store the conductance matrix in a vector
void conduc_mat (const std::vector<int> &nodes, const std::vector<std::vector<std::complex<double>>> &rnode, std::vector<std::complex<double>> &final){
    for (int i = 1; i <= max_node_size(nodes); i++){
        for (int j = 1; j <= max_node_size(nodes); j++){
            std::complex<double> cond = 0;

            std::vector<std::complex<double>> cond_mat;

            intersection_vector(rnode[i-1], rnode[j-1], cond_mat);
            
            for (int v = 0 ; v < cond_mat.size() ; v++){
                cond = cond + cond_mat[v];
            }
            if (i == j){
                final.push_back(cond);
            }
            else{
                final.push_back(-cond);
            }
            
            cond_mat.clear();
        }
    }
}

//Store each nodes of current source into an int vector
std::vector<int> store_nodeI (std::vector<currentSource> node){
    std::vector<int> nodesize;
    for (int i = 0; i < node.size(); i++){
        nodesize.push_back(node[i].n0);
        nodesize.push_back(node[i].n1);
    }
    
    return nodesize;
}

//creates the current vector
std::vector<std::complex<double>> curr_total (std::vector<currentSource> curr, int m){
    std::vector<double> curr_in;
    std::vector<double> curr_out;
    std::vector<std::complex<double>> curr_total;

    for (int i = 0; i <= m; i++){
        curr_total.push_back(0);
        curr_in.push_back(0);
        curr_out.push_back(0);
  
    }

    for (int j = 0; j <= m; j++){
        if (j < curr.size()){
            curr_in[curr[j].n0] += curr[j].dcVal;
            curr_out[curr[j].n1] += -curr[j].dcVal;
        }
    }
 
    for (int n = 0; n <= m; n++){
        curr_total[n] = curr_in[n] + curr_out[n];
    }
    
    return curr_total;
}

//checks for voltage sources
int check_vsource (voltageSource v, int node){
    int ans = 0;
    if (v.n0 == node){
        ans = 1;
    }
    else if (v.n1 == node){
        ans = -1;
    }
    else{
        ans = 0;
    }
    
    return ans;
}

//checks for dependent sources
int check_depV (controlSource g, int col, int row){
    if (g.n0 == row){
        if (g.n2 == col){
            return 1;
        }
        else if (g.n3 == col){
            return (-1);
        }
        else{
            return 0;
        }
    }
    else if (g.n1 == row){
        if (g.n2 == col){
            return -1;
        }
        else if (g.n3 == col){
            return 1;
        }
        else{
            return 0;
        }
    }
    
    return 0;
}

int main(){
    std::vector <std::string> string;
    std::vector <std::string> net;
    std::vector <std::string> netList;
    std::vector <std::complex<double>> output;
    std::vector <std::complex<double>> input;
    std::vector <int> nodes;
    std::string s;
    std::string inputFile;
    
    double points = 0;
    double startFreq = 0;
    double stopFreq = 0;
    double freq = 0;
    int inputNode;
    int outputNode;
    
    std::cout << "Enter name of txt file (include extension): " ;
    std::cin >> inputFile;
    
    std::cout << "Enter input node: ";
    std::cin >> inputNode;
    
    std::cout << "Enter output node: ";
    std::cin >> outputNode;
    
    std::ifstream file (inputFile);
    std::string str;

    if (file.is_open()){
        while (std::getline(file, str)) {
            string.push_back(str);
        }

        std::vector<voltageSource> v;
        std::vector<currentSource> I;
        std::vector<resistor> r;
        std::vector<capacitor> c;
        std::vector<inductor> l;
        std::vector<diode> d;
        std::vector<bjt> q;
        std::vector<mosfet> m;
        std::vector<controlSource> g;
        
        
        for (int i = 0; i < string.size(); i++){
            if (string[i][0] != '.' && string[i][0] != '*' ){
                netList = decodeNet(string[i]);
                store_node(netList, nodes);
                
                if (netList[0][0] == 'V'){
                    voltageSource volt;
                    v.push_back(volt.add(netList));
                }
                else if (netList[0][0] == 'I'){
                    currentSource curr;
                    I.push_back(curr.add(netList));
                }
                else if (netList[0][0] == 'R'){
                    resistor resist;
                    r.push_back(resist.add(netList));
                }
                else if (netList[0][0] == 'C'){
                    capacitor cap;
                    c.push_back(cap.add(netList));
                }
                else if (netList[0][0] == 'L'){
                    inductor in;
                    l.push_back(in.add(netList));
                }
                else if (netList[0][0] == 'D'){
                    diode diode;
                    diode.ssem(diode, netList, r);
                    d.push_back(diode.add(netList));
                }
                else if (netList[0][0] == 'Q'){
                    bjt Bjt;
                    Bjt.ssem(Bjt, g, r, netList);
                    q.push_back(Bjt.add(netList));
                }
                else if (netList[0][0] == 'M'){
                    mosfet mos;
                    m.push_back(mos.add(netList));
                }
                else if (netList[0][0] == 'G'){
                    controlSource cs;
                    g.push_back(cs.add(netList));
                }
            }
            else if (string[i][0] == '.' && string[i][1] == 'a' ) {
                netList = decodeNet(string[i]);
                points = std::stod(netList[2]);
                startFreq = splitVal(netList[3]);
                stopFreq = splitVal(netList[4]);
            }
            else if (string[i] == ".end"){
                break;
            }
        }
        
        freq = 0.000000001;
        double vd = 1;
        double vdnew = 0.7;
        double vbe = 1;
        double vbenew = 1;
        double vce = 1;
        double vcenew = 1;

        do{
            if (d.size() != 0){
                d[0].newton(vd, r);
                vdnew = vd;
            }
            if (q.size() != 0){
                q[0].newtonBJT(vbe, vce, r, g);
                vbenew = vbe;
                vcenew = vce;
            }
            
            std::vector<std::complex<double>> test = overall_imp(r, c, l, freq) ;
            //store which resistors are connected to each node
            std::vector<int> node_list = store_nodercl(r, c, l);
            std::vector<std::vector<std::complex<double>>> rnode;
            Node_vec(node_list, test, rnode, nodes);

            std::vector<std::complex<double>> conmat;
            std::vector<double> conmatreal;
            conduc_mat(nodes, rnode, conmat);
            
            
            for (int o = 0; o <conmat.size(); o++){
                conmatreal.push_back(std::abs(conmat[o]));
            }

            std::vector <std::complex<double>> ci = (curr_total(I, max_node_size(nodes)));
            std::vector <double> currReal;
            
            for (int y = 0; y < ci.size(); y++){
                currReal.push_back(std::abs(ci[y]));
            }

            int k = 0;
            int N = max_node_size(nodes);
            int M = (int)v.size();
            int tot = N + M;

            Eigen::MatrixXd condMatrix(N,N);
            Eigen::VectorXd currVec(N + v.size());
            Eigen::MatrixXd fin(tot,tot);
            Eigen::MatrixXd b(N,M);
            Eigen::MatrixXd C(M,N);
            Eigen::MatrixXd D(M,M);
            
            D << Eigen::MatrixXd::Zero(M,M);

            for (int r = 0; r < v.size(); r++){
                for (int q = 0; q < N; q++){
                    b(q,r) = check_vsource(v[r], q+1);
                }
            }

            C = b.transpose();
            
            if (g.size() != 0){
                for (int row = 0; row < N; row++){
                    currVec(row) = currReal[row+1];
                    for (int col = 0; col < N; col++){
                        condMatrix(row, col) = (conmatreal[k]) + (g[0].transc * check_depV(g[0], col, row));
                        k++;
                    }
                }
            }
            
            else{
                for (int row = 0; row < N; row++){
                    currVec(row) = currReal[row+1];
                    for (int col = 0; col < N; col++){
                        condMatrix(row, col) = (conmatreal[k]);
                        k++;
                    }
                }
            }
            
             for (int a = 0; a < v.size(); a++){
                 if (v[a].ind == 0){
                     currVec(a+N) = v[a].dcVal;
                 }
                 else{
                     currVec(a+N) = std::abs(v[a].complexVal);
                 }
             }
            
            fin << condMatrix, b, C, D;

            Eigen::MatrixXd mat = (fin.inverse()) * currVec;
            
            if (d.size() != 0){
                int node_a = d[0].n0;
                int node_c = d[0].n1;
                vd = mat(node_a) - mat(node_c);
            }
            if (q.size() != 0){
                int node_c = q[0].n0;
                int node_b = q[0].n1;
                int node_e = q[0].n2;
                vbe = mat(node_b) - mat(node_e);
                vce = mat(node_c) - mat(node_e);
            }

            if (d.size() != 0){
                r.pop_back();
            }
            
            if (q.size() != 0){
                r.pop_back();
                r.pop_back();
                g.pop_back();
            }
        }
        while (((int)(vd*10) != (int)(vdnew*10)) && ((int)(vbe*1) != (int)(vbenew*1)) && ((int)(vce*1) != (int)(vcenew*1)));
        
        if (d.size() != 0){
            d[0].newton(vd, r);
        }

        // remake the conductance matrix
        for (int n = 0; freq < stopFreq; n++){
            freq = pow(10,(n/points)) * startFreq;
            std::vector<std::complex<double>> imped = overall_imp(r, c, l, freq) ;
        
            //store which resistors are connected to each node
            std::vector<int> node_list = store_nodercl(r, c, l);
            std::vector<std::vector<std::complex<double>>> rnode;
            Node_vec(node_list, imped, rnode, nodes);

            std::vector<std::complex<double>> conmat ;
            conduc_mat(nodes, rnode, conmat);

            std::vector <std::complex<double>> ci = (curr_total(I, max_node_size(nodes)));

            int k = 0;
            int N = max_node_size(nodes);
            int M = (int)v.size();
            int tot = N + M;

            Eigen::MatrixXcd condMatrix(N,N);
            Eigen::VectorXcd currVec(N + v.size());
            Eigen::MatrixXcd fin(tot,tot);
            Eigen::MatrixXd b(N,M);
            Eigen::MatrixXd C(M,N);
            Eigen::MatrixXd D(M,M);
            D << Eigen::MatrixXd::Zero(M,M);
        
            for (int r = 0; r < v.size(); r++){
                for (int q = 0; q < N; q++){
                    b(q,r) = check_vsource(v[r], q+1);
                }
            }

            C = b.transpose();
            
            if (g.size() != 0){
                for (int row = 0; row < N; row++){
                    currVec(row) = ci[row+1];
                    for (int col = 0; col < N; col++){
                        condMatrix(row, col) = (conmat[k]) + (g[0].transc * check_depV(g[0], col, row));
                        k++;
                    }
                }
            }
            
            else{
                for (int row = 0; row < N; row++){
                    currVec(row) = ci[row+1];
                    for (int col = 0; col < N; col++){
                        condMatrix(row, col) = (conmat[k]);
                        k++;
                    }
                }
            }

             for (int a = 0; a < v.size(); a++){
                 if (v[a].ind == 0){
                     currVec(a+N) = v[a].dcVal;
                 }
                 else{
                     currVec(a+N) = v[a].complexVal;
                 }
             }

            fin << condMatrix, b, C, D;

            Eigen::MatrixXcd mat = (fin.inverse()) * currVec;

            input.push_back(mat(inputNode - 1));
            output.push_back(mat(outputNode - 1));
        }

        std::complex<double> transferFunction;

        for (int count = 0; count < input.size(); count++){
            transferFunction = output[count] / input[count];
            std::cout << transferFunction << std::endl;
        }
        
        file.close();
    }
    
    else{
        std::cout << "Error opening file" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    return 0;
}
