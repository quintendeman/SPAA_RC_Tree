//define generally useful functions, that may be used across project
//ex. bit tricks

#ifndef UTILS_H
#define UTILS_H

static const char PRINT_QUERY = 0;
const double default_epsilon = 0.001;

static bool isNearlyEqual(double a, double b, double epsilon = default_epsilon) {
    return std::abs(a - b) < epsilon;
}


//get index of rightmost 1
template<typename T> 
int r1(T x) {

    return static_cast<int>(log2(x-(x&(x-1))));

}

//get index of leftmost 1
//note: cannot plug bitset in here, must plug in number type
template<typename T>
int l1(T x) {
    if ( x < 0) {
        std::cout << "abort getting leftmost 1 of negative #" << std::endl;
        exit(501);
    }
    return floor(log2(x));
}



/**
 * Extracts a particular bit (counted from the right) from an element
*/
template <typename T>
inline bool extract_bit(T number, int offset_from_right)
{
    return (number >> offset_from_right) & 1;
}

/*
    Returns index of first different bit from the left
    Sets the value of bit in the boolean bit
    Sets the value of b
    sizeof(T) must be less than 16 bytes

    Using bit-tricks, we avoid a for loop
*/
template <typename T>
inline char first_different_bit(const T a, const T b, bool* bit)
{
    T difference = a ^ b;
    if (difference==0) return -1;
    
    char index = l1(difference); //get index of leftmost 1
    if (bit) *bit = extract_bit(b,index); //if bit not a nullptr, store the differing bit in bit
    return index; 

}


//LCA parse input calls
void parse_input(int argc, char* argv[], int& n, int& NUM_TRIALS,int& seed, int& NUM_TREES) {
    for (int i = 1; i < argc; i++) {
        std::string arg=argv[i];
        if (arg=="-n" && i+1 < argc) {
            n=std::stoi(argv[i+1]);

        }
        if (arg=="-trials" && i+1 < argc) {
            NUM_TRIALS=std::stoi(argv[i+1]);
        }
        if (arg=="-seed" && i+1 < argc) {
            seed=std::stoi(argv[i+1]);
        }
        if (arg=="-trees" && i+1 < argc) {
            NUM_TREES=std::stoi(argv[i+1]);
        }

    }
}


void parse_input(int argc, char* argv[], int& n, int& NUM_TRIALS,int& seed, int& pseed, int& NUM_TREES, int& BATCH_SIZE) {
    for (int i = 1; i < argc; i++) {
        std::string arg=argv[i];
        if (arg=="-n" && i+1 < argc) {
            n=std::stoi(argv[i+1]);

        }
        if (arg=="-trials" && i+1 < argc) {
            NUM_TRIALS=std::stoi(argv[i+1]);
        }
        if (arg=="-seed" && i+1 < argc) {
            seed=std::stoi(argv[i+1]);
        }
        if (arg=="-trees" && i+1 < argc) {
            NUM_TREES=std::stoi(argv[i+1]);
        }
        if (arg=="-k" && i+1 < argc) {
            BATCH_SIZE=std::stoi(argv[i+1]);
        }
        if (arg=="-pseed" && i+1 < argc) {
            pseed=std::stoi(argv[i+1]);
        }

    }
}

template<typename T>
void pseq(parlay::sequence<T>& seq) {
    for (int i = 0; i < seq.size(); i++) {
        std::cout << seq[i] << " ";
    }
    std::cout << std::endl;
}


template<typename T>
void pseq(parlay::sequence<T>& seq, std::string message) {
    std::cout << message << ": ";
    for (int i = 0; i < seq.size(); i++) {
        std::cout << seq[i] << " ";
    }
    std::cout << std::endl;
}

#endif