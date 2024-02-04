#pragma warning(disable: 4996)
#ifndef liblgp_hpp
#define liblgp_hpp
#include <numeric>
#include <vector>
#include <set>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <ostream>
#include <functional>
#include <queue>
#include <map>
#include <sstream>
#include <iostream>
#include <math.h>
#include <random>
#include <fstream>
#include <regex>
#ifdef __APPLE__
#include <sys/uio.h>
#include <unistd.h>
#else
#include <io.h>
#include <direct.h>
#include <windows.h>
#include <tchar.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include <locale>
#include <codecvt>
#include <chrono>
#include <thread>

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/matrix_inverse.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/norm.hpp"
#include "glm/gtx/transform.hpp"
#include "glm/gtc/random.hpp"

#include "Eigen/Core"
#include "Eigen/Dense"


//glm modules
//http://glm.g-truc.net/0.9.8/api/modules.html


using namespace std;

namespace liblgp {

    template <typename datum>
    using Vector1 = std::vector<datum>;

    template <typename datum>
    using Vector2 = std::vector<std::vector<datum>>;//two-dimensional array

    template <typename datum>
    using Vector3 = std::vector<std::vector<std::vector<datum>>>;//three-dimensional array

    typedef glm::highp_dvec2 Vector2d;//High precision double precision vector types
    typedef glm::highp_dvec3 Vector3d;
    typedef glm::highp_ivec2 Vector2i;//High precision integer vector types
    typedef glm::highp_ivec3 Vector3i;

    typedef Vector1<Vector2d> Vector2d1;
    typedef Vector2<Vector2d> Vector2d2;
    typedef Vector3<Vector2d> Vector2d3;

    typedef Vector1<Vector3d> Vector3d1;
    typedef Vector2<Vector3d> Vector3d2;
    typedef Vector3<Vector3d> Vector3d3;

    typedef Vector1<bool> Vector1b1;
    typedef Vector2<bool> Vector1b2;
    typedef Vector3<bool> Vector1b3;

    typedef Vector1<int> Vector1i1;
    typedef Vector2<int> Vector1i2;
    typedef Vector3<int> Vector1i3;

    typedef Vector1<double> Vector1d1;
    typedef Vector2<double> Vector1d2;
    typedef Vector3<double> Vector1d3;

    typedef Vector1<std::string> VectorStr1;
    typedef Vector2<std::string> VectorStr2;
    typedef Vector3<std::string> VectorStr3;


    typedef Vector1<Vector2i> Vector2i1;
    typedef Vector2<Vector2i> Vector2i2;
    typedef Vector3<Vector2i> Vector2i3;

    typedef Vector1<Vector3i> Vector3i1;
    typedef Vector2<Vector3i> Vector3i2;
    typedef Vector3<Vector3i> Vector3i3;

    typedef Vector1<std::pair<int, int>> VectorPI1;
    typedef Vector2<std::pair<int, int>> VectorPI2;
    typedef Vector3<std::pair<int, int>> VectorPI3;

    typedef Vector1<std::pair<bool, bool>> VectorPB1;
    typedef Vector2<std::pair<bool, bool>> VectorPB2;
    typedef Vector3<std::pair<bool, bool>> VectorPB3;

    typedef std::tuple<int, int, int> TI3;
    typedef Vector1<std::tuple<int, int, int>> VectorTI3;


    static double DOUBLE_EPSILON = 1.0E-05;
    static double Math_PI = 3.14159265358979323846;
    static double SINGLE_EPSILON = 1.0E-05f;
    static double MAXDOUBLE = 100000000000.0;
    static std::random_device MATHRD;  //Will be used to obtain a seed for the random number engine
    static std::mt19937 MATHGEN(0); //Standard mersenne_twister_engine seeded with rd()
    static std::string CERR_ITER = "  ";

#ifdef __APPLE__
    static bool MACEN = true;
    static bool WINEN = false;
#else
    static bool MACEN = false;
    static bool WINEN = true;
#endif

    struct LiblgpTriMesh
    {
        Vector3d1 vecs;
        Vector1i1 face_id_0;
        Vector1i1 face_id_1;
        Vector1i1 face_id_2;
    };

    struct TimeClock
    {
    public:
        TimeClock() :start(0), end(0), duration(0) { start = clock(); };
        void StartClock() { start = clock(); };
        double EndClock()
        {
            end = clock();
            duration = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;
            return duration;
        };
        double start, end, duration;
    };


    class Functs
    {
    public:

       
#pragma region StatisticsCombinationSet

        template <class Type>
        //Calculate the sample variance after Bessel correction
        static Type Variance(const std::vector<Type>& resultSet)
        {
            Type sum = std::accumulate(std::begin(resultSet), std::end(resultSet), 0.0);
            Type mean = sum / resultSet.size(); 
            Type accum = 0.0;
            std::for_each(std::begin(resultSet), std::end(resultSet), [&](const Type d) {accum += (d - mean) * (d - mean); });
            Type stdev = sqrt(accum / (resultSet.size() - 1)); //Variance

            return stdev;
        }

        //nb>=1
        //Calculate Factorial (Recursive)
        static int Factorial(const int& n)
        {
            if (n > 1)
                return n * Factorial(n - 1);
            else
                return 1;
        };

        //arrangement and combination
        
        static void Combination(const Vector1i3& combs, const int& nb, const Vector1i2& sequence, const int& max_nb, Vector1i3& output)
        {
            if (output.size() > max_nb)return;

            if (nb == combs.size())
            {
                output.emplace_back(sequence);
                return;
            }

            for (int i = 0; i < combs[nb].size(); i++)
            {
                Vector1i2 seq = sequence;
                seq.emplace_back(combs[nb][i]);
                Combination(combs, nb + 1, seq, max_nb, output);
            }
        }

        //arrangement and combination
        //goups={3,2,3}
        //0,0,0 //0,0,1 //0,0,2 //0,1,0 //0,1,1 //0,1,2
        //1,0,0 //1,0,1 //1,0,2 //1,1,0 //1,1,1 //1,1,2
        //2,0,0 //2,0,1 //2,0,2 //2,1,0 //2,1,1 //2,1,2
        static Vector1i2 Selection(const Vector1i1& groups)
        {
            int nb = 1;
            for (auto& group : groups) nb = nb * group;

            Functs::MAssert(nb != 0, "nb == 0 in Selection(const Vector1i1& groups)");

            Vector1i1 nbs(1, nb);
            for (auto& group : groups)
            {
                nbs.emplace_back(nbs.back() / group);
            }
            nbs.erase(nbs.begin());

            Vector1i2 combs;
            for (int i = 0; i < nb; i++)
            {
                int id = i;
                combs.emplace_back(Vector1i1());
                for (auto nbs_ : nbs)
                {
                    int a = (id - id % nbs_) / nbs_;
                    combs.back().emplace_back(a);
                    id = id - a * nbs_;
                }
            }
            return combs;
        }

        //arrangement and combination
        //with repeat selection
        //n=3,m=2
        //1,2
        //1,3
        //2,3
        static Vector1i2 CombNonRepeat(const int& n, const int& m)
        {
            std::vector<std::vector<int>> combs;
            if (n > m)
            {
                vector<int> p, set;
                p.insert(p.end(), m, 1);
                p.insert(p.end(), static_cast<int64_t>(n) - m, 0);
                for (int i = 0; i != p.size(); ++i)
                    set.push_back(i + 1);
                vector<int> vec;
                size_t cnt = 0;
                do {
                    for (int i = 0; i != p.size(); ++i)
                        if (p[i])
                            vec.push_back(set[i]);
                    combs.emplace_back(vec);
                    cnt++;
                    vec.clear();
                } while (prev_permutation(p.begin(), p.end()));
            }
            else
            {
                combs.emplace_back(std::vector<int>());
                for (int i = 1; i <= n; i++)
                    combs.back().emplace_back(i);
            }

            return combs;
        }

        //arrangement and combination
        //N=3,K=2
        //0,0
        //0,1
        //0,2
        //1,1
        //1,2
        //2,2
        static std::vector<std::vector<int>> CombRepeat(const int& N, const int& K)
        {
            auto Combination = [](const int& N, const int& K) {
                std::vector<std::vector<int>> temps;
                std::string bitmask(K, 1); // K leading 1's
                bitmask.resize(N, 0); // N-K trailing 0's
                // print integers and permute bitmask
                do {
                    std::vector<int> temp;
                    for (int i = 0; i < N; ++i) // [0..N-1] integers
                        if (bitmask[i]) temp.emplace_back(i);
                    if (!temp.empty()) temps.emplace_back(temp);
                } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
                return temps;
            };

            std::vector<std::vector<int>> combs = Combination(N + K - 1, K);

            std::vector<std::vector<int>> repeatCombs;
            for (auto& comb : combs) {
                std::vector<int> nbs((int)(static_cast<int64_t>(N) + static_cast<int64_t>(K) - 1), 1);
                for (auto& i : comb) nbs[i] = 0;
                std::vector<int> temp;
                int sum = 0;
                for (auto& nb : nbs) {
                    if (nb == 0)temp.emplace_back(sum);
                    sum += nb;
                }
                if (!temp.empty()) repeatCombs.emplace_back(temp);
            }
            return repeatCombs;
        }

        //remove duplicated elements
        //vec={3,2,1,2,3,4}
        //output:{1,2,3,4}
        static Vector1i1 UniqueSet(const Vector1i1& vec)
        {
            Vector1i1 s;
            for (auto& v : vec)
            {
                if (std::find(s.begin(), s.end(), v) == s.end()) s.emplace_back(v);
            }
            std::sort(s.begin(), s.end());
            return s;
        };
        //Compute Union
        static Vector1i1 SetUnion(const Vector1i1& first, const Vector1i1& second)
        {
            Vector1i1 v = first;
            for (auto& s : second)
                if (std::find(first.begin(), first.end(), s) == first.end())
                    v.emplace_back(s);//Unique to second, join union
            return v;
        }
        //Compress two-dimensional data into one-dimensional and remove duplicates
        static Vector1i1 SetUnion(Vector1i2& sets)
        {
            Vector1i1 start = sets[0];

            for (int i = 1; i < sets.size(); i++)
            {
                start = SetUnion(start, sets[i]);
            }

            return start;
        }
        //Calculate Intersection
        static Vector1i1 SetIntersection(const Vector1i1& first, const Vector1i1& second)
        {
            Vector1i1 v;
            for (auto& s : second)
                if (std::find(first.begin(), first.end(), s) != first.end())
                    v.emplace_back(s);
            return v;
        }
        //The set after subtracting second from first,namely relative complement
        static Vector1i1 SetSubtraction(const Vector1i1& first, const Vector1i1& second)
        {
            Vector1i1 v;
            for (auto& s : first)
                if (std::find(second.begin(), second.end(), s) == second.end())
                    v.emplace_back(s);
            return v;
        }

        //Find a proper subset in a two-dimensional vector so that the elements in the subset contain all the elements in the target, and return the index of that subset
        static Vector1i1 FindSetCombination(Vector1i2& input_)
        {
            auto FSC = [](vector<set<int>>& input, set<int>& target, vector<int>& output)
            {
                set<int> full;
                for (auto it : input) {
                    full.insert(it.begin(), it.end());//Two dimensions reduced to one dimension
                }

                if (!includes(full.begin(), full.end(), target.begin(), target.end())) {
                    return;
                }

                for (int i = static_cast<int>(input.size()) - 1; i > 0; --i) {
                    vector<bool> vec(input.size(), false);
                    fill(vec.begin() + i, vec.end(), true);
                    set<int> comb;

                    do {
                        for (int j = 0; j < vec.size(); ++j) {
                            if (vec[j]) {
                                comb.insert(input[j].begin(), input[j].end());
                            }
                        }

                        if (includes(comb.begin(), comb.end(), target.begin(), target.end())) {
                            for (int j = 0; j < vec.size(); ++j) {
                                if (vec[j]) {
                                    output.push_back(j);
                                }
                            }
                            return;
                        }
                        comb.clear();

                    } while (next_permutation(vec.begin(), vec.end()));
                }
            };


            vector<set<int>> input;
            set<int> target;
            for (auto& a : input_)
            {
                for (int x : a) target.insert(x);
                input.emplace_back(ConvertToSet(a));
            }
            Vector1i1 output;
            FSC(input, target, output);
            return output;
        }



#pragma endregion

        //String data structure
#pragma region StringDataStructure  
        //Integer plus one
        static int INe(const int& i)
        {
            return i + 1;
        }
        //Integer minus one
        static int IPr(const int& i)
        {
            return i - 1;
        }
        //Determine whether a string str contains another string sub
        static bool StringContain(const string& str, const string& sub)
        {
            return str.find(sub) != std::string::npos;
        }
        //Replace String
        //toReplace:The part replaced with
        static std::string StringReplace(const string& source, const string& toReplace, const string& replaceWith)
        {
            size_t pos = 0;
            size_t cursor = 0;
            int repLen = (int)toReplace.length();
            stringstream builder;

            do
            {
                pos = source.find(toReplace, cursor);

                if (string::npos != pos)
                {
                    //copy up to the match, then append the replacement
                    builder << source.substr(cursor, pos - cursor);
                    builder << replaceWith;

                    // skip past the match
                    cursor = pos + repLen;
                }
            } while (string::npos != pos);

            //copy the remainder
            builder << source.substr(cursor);

            return (builder.str());
        }


        //Convert integers to strings
        //p:The number of digits in the converted string, if the original number of digits is insufficient, fill in 0 before it
        static std::string Int2String(int i, int p =0)
        {
            std::stringstream ss;
            std::string str;
            ss << i;
            ss >> str;
            
            if(p>0)
            {
                if(i<0) str = Int2String(abs(i));
               
                auto strc = str.c_str();
                string strp="";
                for(int j=p;j>=1;j--)
                    strp+=j<=str.length()?string(1, strc[j-1]):"0";
                if(i<0)strp="-"+strp;
                return strp;
            }
            
            return str;
        }
        //Convert an array into a string and insert the same string between each element.
        //vecs:The string to convert
        //order:Choose whether to sort
        //insert_str:String to Insert
        template <class Type>
        static std::string Int2String(const std::vector<Type>& vecs, bool order = false, std::string insert_str = "")
        {
            std::vector<Type> a = vecs;
            if (order)sort(a.begin(), a.end());

            std::string str;
            for (int i = 0; i < a.size(); i++)
            {
                str += Int2String(a[i]);
                if (i != a.size() - 1) str += insert_str;
            }
            return str;

        }
        //Convert a two-dimensional array to a string
        template <class Type>
        static std::string Int2String(const std::vector<std::vector<Type>>& vecs, bool order = false, std::string insert_str_0 = "", std::string insert_str_1 = "")
        {
            std::string str;

            for (int i = 0; i < vecs.size(); i++)
            {
                str += IntString(vecs[i], order, insert_str_0);
                if (i != vecs.size() - 1) str += insert_str_1;
            }

            return str;
        }
        //Convert a three-dimensional array to a string
        template <class Type>
        static std::string Int2String(const std::vector<std::vector<std::vector<Type>>>& vecs, bool order = false, std::string insert_str_0 = "", std::string insert_str_1 = "", std::string insert_str_2 = "")
        {
            std::string str;

            for (int i = 0; i < vecs.size(); i++)
            {
                str += IntString(vecs[i], order, insert_str_0, insert_str_1);
                if (i != vecs.size() - 1) str += insert_str_2;
            }

            return str;
        }

        
        //Convert a two-dimensional pair array to a string
        static std::string Int2String(const VectorPI2& vecs, bool order = false,
            const std::string insert_str_0 = "", const std::string insert_str_1 = "", const std::string insert_str_2 = "")
        {
            std::string str;
            for (auto& vec : vecs)
            {
                str += Int2String(vec, order, insert_str_0, insert_str_1) + insert_str_2;
            }
            return str;
        }
        //Convert an array of pair types into a string and insert a string between the two elements of each pair type
        //first+insert_str+second+....
        static std::string Int2String(const VectorPI1& vecs, bool order = false, const std::string insert_str_0 = "", const std::string insert_str_1 = "")
        {
            if (order)
            {
                std::vector<std::pair<int, int>> a = vecs;
                sort(a.begin(), a.end());
                std::string str;
                for (int i = 0; i < a.size(); i++)
                {
                    str += Int2String(a[i].first) + insert_str_0 + Int2String(a[i].second);
                    if (i != a.size() - 1)
                        str += insert_str_1;
                }
                return str;
            }
            else
            {
                std::string str;
                for (int i = 0; i < vecs.size(); i++)
                {
                    str += Int2String(vecs[i].first) + insert_str_0 + Int2String(vecs[i].second);
                    if (i != vecs.size() - 1) str += insert_str_1;
                }
                return str;
            }
        }
        //Convert data of type double to  string
        static std::string Double2String(const double& d, int p = 8)
        {
            double d_ = abs(d);
            double d0 = (int)d_;
            double d1 = d_ - d0;

            {
                double d1_ = d1;
                double d2_ = 0.0;
                if (!Functs::IsAlmostZero(d1_))
                {
                    for (int i = 1; i <= p; i++)
                    {
                        int a = (int)(d1_ * pow(10, i));
                        d1_ = d1_ - (double)a / (double)pow(10, i);
                        d2_ += (double)a / (double)pow(10, i);
                    }

                    auto temp_0 = std::pow(10, -p);
                    auto temp_1 = abs(temp_0 - d1_);

                    if (Functs::IsAlmostZero_Double(temp_1, pow(10, -4)))
                    {
                        d1 = d2_ + temp_0;
                        d_ += (int)d1;
                        d1 = d1 - (int)d1;
                    }
                }
            }

            std::string str;
            str = std::to_string((int)d_) + ".";
            for (int i = 1; i <= p; i++)
            {
                int a = (int)(d1 * pow(10, i));
                str += std::to_string(a);
                d1 = d1 - (double)a / (double)pow(10, i);
            }

            if (d >= 0)
                return str;
            else
                return "-" + str;
        }
        //Convert the double array into a string, where other strings can be inserted between the array elements
        static std::string Double2String(const Vector1d1& ds_, int p = 8, bool order = false, const std::string insert_str_0 = "")
        {
            Vector1d1 ds = ds_;
            if (order)std::sort(ds.begin(), ds.end());
            std::string str;
            for (int i = 0; i < ds.size(); i++)
                str += Double2String(ds[i], p) + (i == ds.size() - 1 ? "" : insert_str_0);
            return str;
        }
        //Convert an array to a string
        template <class Type>
        static std::string Vector2String(const Type& v, const string insert_str = "", int p = 3)
        {
            std::string str;
            for (int i = 0; i < v.length(); i++)
            {
                if (i == v.length() - 1)
                    str += Double2String(v[i], p);
                else
                    str += Double2String(v[i], p) + insert_str;
            }
            return str;
        };
        //Sort a set of three-dimensional vectors and convert them into a string(array dimension:one-dimensional)
        static std::string Vector2String(const Vector3d1& vecs)
        {
            auto Comp = [](Vector3d& v_0, Vector3d& v_1)
            {
                if (IsAlmostZero(abs(v_1[0] - v_0[0])))
                {
                    if (IsAlmostZero(abs(v_1[1] - v_0[1])))
                    {
                        if (IsAlmostZero(abs(v_1[2] - v_0[2])))
                            return true;
                        return v_0[2] < v_1[2];
                    }
                    return v_0[1] < v_1[1];
                }
                return v_0[0] < v_1[0];
            };

            Vector3d1 vecs_1 = vecs;
            std::sort(vecs_1.begin(), vecs_1.end(), Comp);

            std::string str;
            for (auto& p : vecs_1)
            {
                double x = floor(p[0] * 10.0f + 0.5) / 10.0f;
                double y = floor(p[1] * 10.0f + 0.5) / 10.0f;
                double z = floor(p[2] * 10.0f + 0.5) / 10.0f;
                str += Double2String(x);
                str += Double2String(y);
                str += Double2String(z);
            }

            return str;
        }


        //Sort a set of three-dimensional vectors and convert them into a string(array dimension:three-dimensional)
        static std::string Vector2String(const Vector3d3& vecs_3)
        {
            Vector3d1 vecs_1;
            for (int i = 0; i < vecs_3.size(); i++)
                for (int j = 0; j < vecs_3[i].size(); j++)
                    for (int k = 0; k < vecs_3[i][j].size(); k++)
                        vecs_1.emplace_back(vecs_3[i][j][k]);

            return Vector2String(vecs_1);
        };
        //Convert string to double  data
        static double String2Double(const string& str)
        {
            istringstream iss(str);
            double num;
            iss >> num;
            return num;
        }
        //Convert string to type（custom type or basic data type） data
        template <class Type>
        static Type String2Num(const string& str)
        {
            istringstream iss(str);
            Type num;
            iss >> num;
            return num;
        }
        //Encode a null-terminated character array,convert a string to a unique unsigned integer encoding (C style)
        static constexpr unsigned int StringEncode(const char* str, int h = 0)
        {
            return !str[h] ? 5381 : (StringEncode(str, h + 1) * 33) ^ str[h];
        }
      
        static constexpr unsigned int StringEncode(const string& str, int h = 0)
        {
            return StringEncode(str.c_str(), h);
        }
        //Splitting the string str by  delimiter(delimiter is  char)
        static vector<string> SplitStr(const string& str, const char& delimiter)
        {
            vector<string> internal_strs;
            stringstream ss(str); // Turn the string into a stream.
            string tok;
            while (getline(ss, tok, delimiter))
                internal_strs.emplace_back(tok.c_str());
            return internal_strs;
        }
        //Splitting the string str by  delimiter(delimiter is string)
        static vector<string> SplitStr(const string& str, const string& delimiter)
        {
            vector<string> internal_strs;
            std::string s = str;
            size_t pos = 0;
            std::string token;
            while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                internal_strs.push_back(token);
                s.erase(0, pos + delimiter.length());
            }
            internal_strs.push_back(s);
            return internal_strs;
        }
        //Separate strings by delimiter and convert the separated strings into double data(delimiter is  char)
        static vector<double> SplitD(const string& str, const char& delimiter)
        {
            vector<double> internal_d;
            stringstream ss(str); // Turn the string into a stream.
            string tok;
            while (getline(ss, tok, delimiter))
                internal_d.emplace_back(atof(tok.c_str()));
            return internal_d;
        };
        //Separate strings by delimiter and convert the separated strings into integer data(delimiter is  char)
        static vector<int> SplitI(const string& str, const char& delimiter) {
            vector<int> internal_d;
            stringstream ss(str); // Turn the string into a stream.
            string tok;
            while (getline(ss, tok, delimiter))
                internal_d.emplace_back(atoi(tok.c_str()));
            return internal_d;
        };
        //Return a specified size integer array and shuffle the order of its elements
        //size::the size of array
        static Vector1i1 ShuffleVector(const int& size)
        {
            Vector1i1 shuffle_vec;
            for (int i = 0; i < size; i++)
                shuffle_vec.emplace_back(i);
            std::shuffle(shuffle_vec.begin(), shuffle_vec.end(), MATHGEN);
            return shuffle_vec;
        };
        //Return an increasing array of integers with a range of [minI, maxI]
        static Vector1i1 IncreaseVector(const int& minI, const int& maxI)
        {
            MAssert(minI >= 0 && maxI >= 0 && minI <= maxI, "minI>=0&&maxI>=0&&minI<=maxI");
            Vector1i1 vec;
            for (int i = minI; i <= maxI; i++) vec.push_back(i);
            return vec;
        }
        //Return a decreasing array of integers with a range of [minI, maxI]
        static Vector1i1 DecreaseVector(const int& maxI, const int& minI)
        {
            MAssert(minI >= 0 && maxI >= 0 && minI <= maxI, "minI>=0&&maxI>=0&&minI<=maxI");
            Vector1i1 vec;
            for (int i = minI; i <= maxI; i++) vec.push_back(i);
            std::reverse(vec.begin(), vec.end());
            return vec;
        }
        //Convert a vector type to a set type
        static set<int> ConvertToSet(const vector<int>& v)
        {
            set<int> s;
            for (int x : v) s.insert(x);
            return s;
        };
        //Convert a set type to a vector type
        static vector<int> ConvertToVector(const set<int>& s)
        {
            vector<int> v;
            for (auto x : s)v.emplace_back(x);
            return v;
        }
        //Find value based on key
        template <class T1, class T2>
        static T2 MapFind(const std::map<T1, T2>& mt, const T1& t1)
        {
            if (mt.find(t1) == mt.end())
                Functs::MAssert("The input is not in this map.");
            return mt.find(t1)->second;
        }
        //Find if there is a value corresponding to the key in the map data structure
        template <class T1, class T2>
        static bool MapContain(const std::map<T1, T2>& mt, const T1& t1)
        {
            return mt.find(t1) != mt.end();
        }
        //Remove duplicate elements from integer arrays and sort them in ascending order
        static Vector1i1 RemoveDuplicate(const Vector1i1& vec_)
        {
            Vector1i1 vec = vec_;
            sort(vec.begin(), vec.end());
            vec.erase(unique(vec.begin(), vec.end()), vec.end());//Unique will move the repeating element to the end and return the iterator for the first repeating element
            return vec;
        }
        //Convert Vector3d type in the Eigen library to custom Vector3d type 
        static Vector3d Eigen2Vector(const Eigen::Vector3d& vec)
        {
            return Vector3d(vec[0], vec[1], vec[2]);
        };
        //Convert an array with elements of Eigen:: Vector3d type to a custom Vector3d1 type  
        static Vector3d1 Eigen2Vector(const std::vector<Eigen::Vector3d>& vecs)
        {
            Vector3d1 vs;
            for (auto& vec : vecs)
                vs.push_back(Eigen2Vector(vec));
            return vs;
        };

#pragma endregion

#pragma region BasicGeomFunctions
        //Get the length of a vector
        template <class Type>
        static double GetLength(const Type& v) {
            return glm::length(v);
        }
        //Calculate the angle between two vectors
        template <class Type>
        static double GetAngleBetween(const Type& v1, const Type& v2) {
            double d = glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2));//点乘计算夹角余弦
            if (IsAlmostZero(d - 1.0))
                return 0.0;
            if (IsAlmostZero(d + 1.0))
                return Math_PI;
            return glm::acos(d);
        }
        //Convert radians to angles
       
        static double Radian2Angle(const double& radian)
        {
            return radian / Math_PI * 180.0;
        }
        //Convert angles to radians
        static double Angle2Radian(const double& angle)
        {
            return angle / 180.0 * Math_PI;
        }

        //Scaling a 3D vector, scale the length of vector v to the specified length 
        static Vector3d GetVectorLength(const Vector3d& v, const double& length)
        {
            if (IsAlmostZero(length)) return Vector3d(0.0, 0.0, 0.0);//Both scaling factor and vector length cannot be zero
            double l = GetLength(v);
            MAssert(!IsAlmostZero(l), "SetVectorLength: input length is zero.");
            return Vector3d(v[0] / l * length, v[1] / l * length, v[2] / l * length);
        }
        //Scaling a 2D vector
        static Vector2d GetVectorLength(const Vector2d& v, const double& length)
        {
            if (IsAlmostZero(length)) return Vector3d(0.0, 0.0, 0.0);
            double l = GetLength(v);
            MAssert(!IsAlmostZero(l), "SetVectorLength: input length is zero.");
            return Vector2d(v[0] / l * length, v[1] / l * length);
        }
        //Set the length of a 3D vector
        static void SetVectorLength(Vector3d& v, const double& length)
        {
            v = GetVectorLength(v, length);
        }
        //Set the length of a 2D vector
        static void SetVectorLength(Vector2d& v, const double& length)
        {
            v = GetVectorLength(v, length);
        }



        //existing bugs in this function
        //Finding a vector perpendicular to vector v
        static Vector3d Vector3dBase(const Vector3d& v)
        {
            Vector3d n(1.0, 1.0, 1.0);
            if (!IsAlmostZero(v[0])) {
                n[0] = -(v[1] + v[2]) / v[0];//v[0]*n[0]+v[1]*1+v[2]*1=0
                return n;
            }
            if (!IsAlmostZero(v[1])) {
                n[1] = -(v[0] + v[2]) / v[1];
                return n;
            }
            if (!IsAlmostZero(v[2])) {
                n[2] = -(v[0] + v[1]) / v[2];
                return n;
            }
            return n;
        }
        //Find the normal vector of a plane determined by v1,v2
        static Vector3d GetNormal(const Vector3d& v0, const Vector3d& v1)
        {
            MAssert(!IsAlmostZero(GetLength(v0)), "v0 is zero.");
            MAssert(!IsAlmostZero(GetLength(v1)), "v1 is zero.");
            double angle = GetAngleBetween(v0, v1);
            if (IsAlmostZero(angle) || IsAlmostZero(angle - Math_PI))return Vector3dBase(v0);////Two vectors are collinear, finding any perpendicular vector is the normal vector
            return GetCrossproduct(v0, v1);//If not collinear, cross multiply to find the normal vector
        }
        //Find the center point of a set of three-dimensional points(array dimension: one-dimensional)
        static Vector3d GetCenter(const Vector3d1& points)
        {
            Vector3d center(0.0, 0.0, 0.0);
            for (int i = 0; i < points.size(); i++)
                center += points[i];
            center = center / (double)points.size();
            return center;
        }
        //Find the center point of a set of three-dimensional points(array dimension: two-dimensional)
        static Vector3d GetCenter(const Vector3d2& points)
        {
            Vector3d center(0.0, 0.0, 0.0);
            int nb = 0;
            for (int i = 0; i < points.size(); i++)
            {
                for (int j = 0; j < points[i].size(); j++)
                    center += points[i][j];
                nb += static_cast<int>(points[i].size());
            }
            center = center / (double)nb;
            return center;
        }
        //Find the center point of a set of three-dimensional points(array dimension: three-dimensional)
        static Vector3d GetCenter(const Vector3d3& points)
        {
            Vector3d center(0.0, 0.0, 0.0);
            int nb = 0;
            for (int i = 0; i < points.size(); i++)
            {
                for (int j = 0; j < points[i].size(); j++)
                {
                    for (int k = 0; k < points[i][j].size(); k++)
                    {
                        center += points[i][j][k];
                        nb++;
                    }
                }
            }
            center = center / (double)nb;
            return center;
        }
        //Find the center point of a set of two-dimensional points(array dimension: one-dimensional)
        static Vector2d GetCenter(const std::vector<Vector2d>& points)
        {
            Vector2d center(0.0, 0.0);
            for (int i = 0; i < points.size(); i++)
                center += points[i];
            center = center / (double)points.size();
            return center;
        }
        //Find the center point of a set of two-dimensional points(array dimension: two-dimensional)
        static Vector2d GetCenter(const std::vector<std::vector<Vector2d>>& points)
        {
            Vector2d center(0.0, 0.0);
            int nb = 0;
            for (int i = 0; i < points.size(); i++)
            {
                for (int j = 0; j < points[i].size(); j++)
                    center += points[i][j];
                nb += static_cast<int>(points[i].size());
            }
            center = center / (double)nb;
            return center;
        }
        //Find the bounding box of a set of three-dimensional points(array dimension: one-dimensional)
        static void GetBoundingBox(const Vector3d1& points, Vector3d& minimal_corner, Vector3d& maximal_corner)
        {
            minimal_corner = Vector3d(MAXDOUBLE, MAXDOUBLE, MAXDOUBLE);
            maximal_corner = Vector3d(-MAXDOUBLE, -MAXDOUBLE, -MAXDOUBLE);

            for (int i = 0; i < points.size(); i++)
            {
                minimal_corner[0] = min(minimal_corner[0], points[i][0]);
                minimal_corner[1] = min(minimal_corner[1], points[i][1]);
                minimal_corner[2] = min(minimal_corner[2], points[i][2]);
                maximal_corner[0] = max(maximal_corner[0], points[i][0]);
                maximal_corner[1] = max(maximal_corner[1], points[i][1]);
                maximal_corner[2] = max(maximal_corner[2], points[i][2]);
            }
        }

        //Find the bounding box of a set of three-dimensional points(array dimension: two-dimensional)
        static void GetBoundingBox(const Vector3d2& points, Vector3d& minimal_corner, Vector3d& maximal_corner)
        {
            //Two points can determine the range of the bounding box
            minimal_corner = Vector3d(MAXDOUBLE, MAXDOUBLE, MAXDOUBLE);//Bottom left corner point of bounding box (minimum coordinate point)
            maximal_corner = Vector3d(-MAXDOUBLE, -MAXDOUBLE, -MAXDOUBLE);//Upper right corner point of bounding box (maximum coordinate point)
            for (int i = 0; i < points.size(); i++)
            {
                for (int j = 0; j < points[i].size(); j++)
                {
                    minimal_corner[0] = min(minimal_corner[0], points[i][j][0]);
                    minimal_corner[1] = min(minimal_corner[1], points[i][j][1]);
                    minimal_corner[2] = min(minimal_corner[2], points[i][j][2]);
                    maximal_corner[0] = max(maximal_corner[0], points[i][j][0]);
                    maximal_corner[1] = max(maximal_corner[1], points[i][j][1]);
                    maximal_corner[2] = max(maximal_corner[2], points[i][j][2]);
                }
            }
        }
       
        //Find the bounding box of a set of three-dimensional points(array dimension: one-dimensional)
        static void GetBoundingBox(const std::vector<Vector2d>& points, Vector2d& minimal_corner, Vector2d& maximal_corner)
        {
            minimal_corner = Vector2d(MAXDOUBLE, MAXDOUBLE);
            maximal_corner = Vector2d(-MAXDOUBLE, -MAXDOUBLE);
            for (int i = 0; i < points.size(); i++)
            {
                minimal_corner[0] = min(minimal_corner[0], points[i][0]);
                minimal_corner[1] = min(minimal_corner[1], points[i][1]);
                maximal_corner[0] = max(maximal_corner[0], points[i][0]);
                maximal_corner[1] = max(maximal_corner[1], points[i][1]);
            }
        }
        //Find the bounding box of a set of three-dimensional points(array dimension: two-dimensional)
        static void GetBoundingBox(const std::vector<std::vector<Vector2d>>& points, Vector2d& minimal_corner, Vector2d& maximal_corner)
        {
            minimal_corner = Vector2d(MAXDOUBLE, MAXDOUBLE);
            maximal_corner = Vector2d(-MAXDOUBLE, -MAXDOUBLE);
            for (int i = 0; i < points.size(); i++)
            {
                for (int j = 0; j < points[i].size(); j++)
                {
                    minimal_corner[0] = min(minimal_corner[0], points[i][j][0]);
                    minimal_corner[1] = min(minimal_corner[1], points[i][j][1]);
                    maximal_corner[0] = max(maximal_corner[0], points[i][j][0]);
                    maximal_corner[1] = max(maximal_corner[1], points[i][j][1]);
                }
            }
        }
        //Calculate the radius of the outer circle of a triangle
        static double CircumCircleRaidius(const Vector2d& v0, const Vector2d& v1, const Vector2d& v2)
        {
            double a = GetLength(v0 - v1);//edge length
            double b = GetLength(v0 - v2);
            double c = GetLength(v1 - v2);
            double p = (a + b + c) / 2.0;
            double area = (4.0 * pow(p * (p - a) * (p - b) * (p - c), 0.5));//Helen's formula for calculating area
            double radius;

            if (IsAlmostZero(area))
            {
                double max_l = a;
                max_l = max(max_l, b);
                max_l = max(max_l, c);
                radius = 10 * max_l;//The area is approximately 0, which can be considered as a very linear triangle. Set the area to 10 times the longest side
            }
            else
                radius = a * b * c / area;//Calculating the radius of a circumscribed circle from the area
            return radius;
        }
        //Calculate the area of a triangle
        static double GetTriangleArea(const Vector3d& v0, const Vector3d& v1, const Vector3d& v2)
        {
            double a = GetDistance(v0, v1);
            double b = GetDistance(v2, v1);
            double c = GetDistance(v0, v2);
            double p = (a + b + c) / 2.0;
            return sqrt(p * (p - a) * (p - b) * (p - c));//Heron's formula
        }
        //Distance between two points
        template <class Type>
        static double GetLength(const Type& v0, const Type& v1) {
            return GetLength(v0 - v1);
        }

        template <class Type>
        static double GetDistance(const Type& v0, const Type& v1) {
            return GetLength(v0 - v1);
        }
        //Distance from closest point
        template <class Type>
        static double GetDistance(const Type& v, const std::vector<Type>& vs)
        {
            double min_d = MAXDOUBLE;
            for (const auto& iter : vs)
                min_d = min(min_d, GetDistance(v, iter));
            return min_d;
        }
        
        template <class Type>
        static double GetDistance(const Type& v, const std::vector<std::vector<Type>>& vs)
        {
            double min_d = MAXDOUBLE;
            for (const auto& iter : vs)
                min_d = min(min_d, GetDistance(v, iter));
            return min_d;
        }
        //The average distance between two sets of points
        template <class Type>
        static double GetDistance(const std::vector<std::vector<Type>>& vs1, const std::vector<std::vector<Type>>& vs2)
        {
            double total_d = 0.0;
            int nb = 0;
            for (auto& vec : vs1)
            {
                for (auto& v : vec)
                {
                    total_d = total_d + GetDistance(v, vs2);
                    nb++;
                }
            }
            if (nb != 0) total_d = total_d / (double)nb;
            return total_d;
        }

        template <class Type>
        static double GetDistanceNONE(
            const std::vector<std::vector<Type>>& vs1,
            const std::vector<std::vector<Type>>& vs2)
        {
            return (GetDistance(vs1, vs2) + GetDistance(vs2, vs1)) / 2.0;
        }
        //Find the index of the closest point
        template <class Type>
        static int GetNearestPointIndex(const Type& v, const std::vector<Type>& vs)
        {
            int index = -1;
            double min_d = MAXDOUBLE;
            for (int i = 0; i < vs.size(); i++)
            {
                double cur_d = GetDistance(v, vs[i]);
                if (cur_d < min_d)
                {
                    min_d = cur_d;
                    index = i;
                }
            }
            return index;
        }
        //Calculate the distance between adjacent vectors in a set of vectors and sum it (2D vector)
        static double GetLength(const std::vector<Vector2d>& points)
        {
            double length = 0.0;

            for (int i = 0; i < points.size(); i++)
                length += GetLength(points[i], points[(static_cast<int64_t>(i) + 1) % points.size()]);

            return length;
        }
        //Calculate the distance between adjacent vectors in a set of vectors and sum it (3D vector)
        static double GetLength(const Vector3d1& points)
        {
            double length = 0.0;

            for (int i = 0; i < points.size(); i++)
                length += GetLength(points[i], points[(static_cast<int64_t>(i) + 1) % points.size()]);

            return length;
        }
        //Remove points with a distance less than 0.00001 between adjacent points
        static Vector2d1 RemoveClosePoints(const Vector2d1& xys_)
        {
            Vector2d1 xys = xys_;
            std::vector<int> remove_int;//Store the index of the vector to be deleted
            if (xys.size() > 2)
            {
                for (int i = 0; i < xys.size() - 1; i++)
                {
                    double d = GetDistance(xys[i], xys[(int)(i + 1)]);//Calculate the distance between adjacent vectors
                    if (d < 0.00001) remove_int.push_back(i + 1);
                }

                for (int i = (int)(remove_int.size() - 1); i >= 0; i--)
                {
                    xys.erase(xys.begin() + remove_int[i]);
                }
            }
            return xys;
        }
        //Project a point in 3D space onto a plane.Return projection point coordinates
        //p:Points to project
        //planar_location:A point in the plane
        //planar_direction:A normal vector of a plane
        static Vector3d PlaneProject(const Vector3d& planar_location, const Vector3d& planar_direction, const Vector3d& p)
        {
            if (IsAlmostZero(GetLength(planar_location, p)))
                return planar_location;//In this case, the projection point is planar_ location

            double angle = GetAngleBetween(planar_direction, p - planar_location);//The angle between the normal vector and the vector formed by the point to be projected and a known point in the plane
            double length = GetLength(planar_location, p);//The distance between a known point and the point to be projected

            auto backup_planar_direction = planar_direction;

            if (angle <= Math_PI / 2.0)
                return p - GetVectorLength(backup_planar_direction, length * sin(Math_PI / 2.0 - angle));//Move the point to be projected towards the normal direction to obtain the projection coordinates
            else
                return p + GetVectorLength(backup_planar_direction, length * sin(angle - Math_PI / 2.0));

        }
        //Generate a three-dimensional vector with a random direction
        static Vector3d GetRandomDirection(const double& direction_length = 1.0)
        {
            double alpha_angle = rand() / double(RAND_MAX) * 2.0 * liblgp::Math_PI;
            double alpha_beta = rand() / double(RAND_MAX) * 2.0 * liblgp::Math_PI;
            auto direction_0 = liblgp::Functs::RotationAxis(Vector3d(direction_length, 0.0, 0.0), alpha_angle, Vector3d(0.0, 1.0, 0.0));
            auto direction_axis = liblgp::Functs::GetCrossproduct(direction_0, Vector3d(0.0, 1.0, 0.0));
            return liblgp::Functs::RotationAxis(direction_0, alpha_beta, direction_axis);
        }



        //https://medium.com/@all2one/generating-uniformly-distributed-points-on-sphere-1f7125978c4c

        //method: NormalDeviate; TrigDeviate; CoordinateApproach;MinimalDistance;RegularDistribution;
        //Generate a set of random vectors with different distributions based on different commands
        static Vector3d1 GetRandomDirections(const int& count, const string& method)
        {
            if (method == "NormalDeviate") return GetRandomDirections_Normal_Deviate(count);
            if (method == "TrigDeviate") return GetRandomDirections_Trig_method(count);
            if (method == "CoordinateApproach") return GetRandomDirections_Coordinate_Approach(count);
            if (method == "MinimalDistance") return GetRandomDirections_Minimal_Distance(count);
            if (method == "RegularDistribution") return GetRandomDirections_Regular_Distribution(count);

            MAssert("Input method does not be implemented: " + method);
            return Vector3d1();
        }
        //Using the method of normal distribution to generate random unit vectors with regular distributions
        static Vector3d1 GetRandomDirections_Regular_Distribution(const int& count)
        {
            Vector3d1 directions;
            double a = 4.0 * Math_PI * 1.0 / static_cast<double>(count);
            double d = sqrt(a);
            size_t num_phi = (size_t)round(Math_PI / d);
            double d_phi = Math_PI / static_cast<double>(num_phi);
            double d_theta = a / d_phi;
            for (int m = 0; m < num_phi; ++m) {
                double phi = Math_PI * (m + 0.5) / num_phi;
                size_t num_theta = (size_t)round(2 * Math_PI * sin(phi) / d_theta);
                for (int n = 0; n < num_theta; ++n) {
                    double theta = 2 * Math_PI * n / static_cast<double>(num_theta);
                    Vector3d p;
                    p.x = sin(phi) * cos(theta);
                    p.y = sin(phi) * sin(theta);
                    p.z = cos(phi);
                    directions.push_back(p);
                }
            }
            return directions;
        }

        //Generating random unit vectors using normal distribution
        static Vector3d1 GetRandomDirections_Normal_Deviate(const int& count)
        {
            std::mt19937 rnd;
            std::normal_distribution<double> dist(0.0, 1.0);//Random numbers with normal distribution, range 0-1

            Vector3d1 directions;
            for (int i = 0; i < count; ++i)
            {
                bool bad_luck = false;
                do
                {
                    double x = dist(rnd);
                    double y = dist(rnd);
                    double z = dist(rnd);
                    double r2 = x * x + y * y + z * z;
                    if (r2 == 0)
                        bad_luck = true;
                    else
                    {
                        bad_luck = false;
                        double r = sqrt(r2);
                        directions.push_back(Vector3d(x / r, y / r, z / r));
                    }
                } while (bad_luck);
            }

            return directions;
        }
        //Using trigonometric function method to generate random unit vectors
        static Vector3d1 GetRandomDirections_Trig_method(const int& count)
        {
            std::mt19937 rnd;
            std::uniform_real_distribution<double> dist(0.0, 1.0);

            Vector3d1 directions;
            for (int i = 0; i < count; ++i)
            {
                double z = 2.0 * dist(rnd) - 1.0;//The range of z is -1 to 1
                double t = 2.0 * Math_PI * dist(rnd);//The range of t is 0 to 2 π
                double r = sqrt(1.0 - z * z);
                directions.push_back(Vector3d(r * cos(t), r * sin(t), z));
            }
            return directions;
        };
        //Randomly sampling a set of vectors within a unit sphere
        static Vector3d1 GetRandomDirections_Coordinate_Approach(const int& count)
        {
            std::mt19937 rnd;
            std::uniform_real_distribution<double> dist(-1.0, 1.0);//The range of random numbers is [-1,1]

            Vector3d1 directions;
            for (int i = 0; i < count; ++i)
            {
                bool rejected = false;//Determine if the generated vector is within the unit sphere
                do
                {
                    double u = dist(rnd);
                    double v = dist(rnd);
                    double s = u * u + v * v;
                    if (s > 1.0)
                        rejected = true;
                    else
                    {
                        rejected = false;
                        double a = 2.0 * sqrt(1.0 - s);
                        directions.push_back(Vector3d(a * u, a * v, 2.0 * s - 1.0));
                    }
                } while (rejected);
            }
            return directions;
        }

        //random sample a set of directions on the Gaussian Sphere
        
        static Vector3d1 GetRandomDirections_Minimal_Distance(const int& dns, const int dis_iters = 100)
        {
            double gaussion_sphere_radius = 1.0;
            double idea_distance = 2 * gaussion_sphere_radius / sqrt(dns);//Ideal distance, determined based on Gaussian sphere radius and dns

            Vector3d1 directions;
            for (int i = 0; i < dns; i++)
            {
                OutputIterInfo("Random Directions", dns, i, 10);

                for (int j = 0; j < dis_iters; j++)
                {
                    //glm::ballRand(gaussion_sphere_radius) is slower than my solution
                    Vector3d random_direction = GetRandomDirection(gaussion_sphere_radius);

                    auto dis = Functs::GetDistance(random_direction, directions);//Calculate the minimum distance between the current randomly generated vector and the generated vector
                    if (dis > idea_distance)//If the distance is large enough, add it to the generated vector sequence
                    {
                        directions.push_back(random_direction);
                        break;
                    }
                }

            }
            return directions;
        }
        //Group the segments  according to their connection relationships, and store the segments in each group in lines
        static void ConnectingSegments(const Vector3d2& segments, Vector3d2& lines)
        {
            //save connecting relations
            std::vector<bool> used(segments.size(), false);
            std::vector<int> relations;
#pragma region get_relations
            for (int i = 0; i < segments.size(); i++)
            {
                for (int j = i + 1; j < segments.size(); j++)
                {
                    if (i != j && !used[i] && !used[j])
                    {
                        double l_0_0 = GetLength(segments[i][0], segments[j][0]);
                        double l_0_1 = GetLength(segments[i][0], segments[j][1]);
                        double l_1_0 = GetLength(segments[i][1], segments[j][0]);
                        double l_1_1 = GetLength(segments[i][1], segments[j][1]);

                        bool b_0_0 = IsAlmostZero_Double(l_0_0, DOUBLE_EPSILON);
                        bool b_0_1 = IsAlmostZero_Double(l_0_1, DOUBLE_EPSILON);
                        bool b_1_0 = IsAlmostZero_Double(l_1_0, DOUBLE_EPSILON);
                        bool b_1_1 = IsAlmostZero_Double(l_1_1, DOUBLE_EPSILON);

                        if ((b_0_0 && b_1_1) || (b_0_1 && b_1_0))
                        {
                            used[j] = true;
                            continue;
                        }

                        if (b_0_0)
                        {
                            relations.push_back(i);
                            relations.push_back(0);
                            relations.push_back(j);
                            relations.push_back(0);
                            continue;
                        }
                        if (b_0_1)
                        {
                            relations.push_back(i);
                            relations.push_back(0);
                            relations.push_back(j);
                            relations.push_back(1);
                            continue;
                        }
                        if (b_1_0)
                        {
                            relations.push_back(i);
                            relations.push_back(1);
                            relations.push_back(j);
                            relations.push_back(0);
                            continue;
                        }
                        if (b_1_1)
                        {
                            relations.push_back(i);
                            relations.push_back(1);
                            relations.push_back(j);
                            relations.push_back(1);
                            continue;
                        }
                    }
                }
            }
#pragma endregion

            std::vector<std::vector<int>> ones;


            while (true)
            {
                int index = -1;
                int end = -1;

                for (int i = 0; i < segments.size(); i++)
                {
                    if (!used[i]) {
                        index = i;
                        end = 0;
                        used[i] = true;
                        break;
                    }
                }

                if (index < 0)break;

                Vector3d1 line(1, segments[index][end]);

                std::vector<int> one(1, index);

                while (true)
                {
                    end = 1 - end;
                    bool search = false;
                    for (int i = 0; i < relations.size(); i = i + 4)
                    {
                        if (relations[i] == index && relations[static_cast<int64_t>(i) + 1] == end && !used[relations[static_cast<int64_t>(i) + 2]])
                        {
                            line.push_back(segments[relations[static_cast<int64_t>(i) + 2]][relations[static_cast<int64_t>(i) + 3]]);
                            one.push_back(relations[static_cast<int64_t>(i) + 2]);
                            index = relations[static_cast<int64_t>(i) + 2];
                            end = relations[static_cast<int64_t>(i) + 3];
                            used[index] = true;
                            search = true;
                            break;
                        }
                        if (relations[static_cast<int64_t>(i) + 2] == index && relations[static_cast<int64_t>(i) + 3] == end && !used[relations[i]])
                        {
                            line.push_back(segments[relations[i]][relations[static_cast<int64_t>(i) + 1]]);
                            one.push_back(relations[i]);
                            index = relations[i];
                            end = relations[static_cast<int64_t>(i) + 1];
                            used[index] = true;
                            search = true;
                            break;
                        }
                    }
                    if (!search) { break; }
                }

                ones.push_back(one);
                lines.push_back(line);
            }
        }
        //Calculate the intersection point between a plane and a ray
        //planar_location:Known point in the plane
        //planar_direction：Plane normal 
        //vectorray_location：Ray origin
        //ray_vector：Ray Direction
        static Vector3d IntersectPointPlane2Ray(const Vector3d& planar_location, Vector3d& planar_direction,
            const Vector3d& ray_location, Vector3d& ray_vector)
        {
            Vector3d project_point = PlaneProject(planar_location, planar_direction, ray_location);//The projection point of the starting point of a ray on a plane
            double distance = GetDistance(ray_location, project_point);//The distance between the starting point of the ray and the projection point
            if (IsAlmostZero(GetLength(project_point, ray_location)))
                return ray_location;
            double angle = GetAngleBetween(ray_vector, project_point - ray_location);//The angle between the direction of a ray and the vector formed by the starting point and projection point of the ray
            double length = distance / cos(angle);//Distance between starting point and intersection point
            return ray_location + GetVectorLength(ray_vector, length);//Starting point+vector direction * length
        }
        //Calculate the unit normal vector of a polygon face
        static Vector3d ComputeNormalFromPolyline(const Vector3d1& points)
        {
            Vector3d planar_direction;
            planar_direction = GetNormal(points[0] - points[1], points[2] - points[1]);//normal vector
            SetVectorLength(planar_direction, 1.0);//normalization
            return planar_direction;
        }
        //Calculate the unit normal vector and the position of a polygonal plane 
        static void  ComputePlanarFromPolyline(Vector3d& planar_location, Vector3d& planar_direction, const Vector3d1& points)
        {
            planar_location = points[0];//The coordinates of a point in a plane
            planar_direction = GetNormal(points[0] - points[1], points[2] - points[1]);
            SetVectorLength(planar_direction, 1.0);//unit normal vector
        }



        // Compute barycentric coordinates (u, v, w) for point p with respect to triangle (a, b, c)
        
        static void Barycentric(const Vector3d& p, const Vector3d& a, const Vector3d& b, const Vector3d& c, double& u, double& v, double& w)
        {
            Vector3d v0 = GetMinus(b, a), v1 = GetMinus(c, a), v2 = GetMinus(p, a);

            double d00 = GetDotproduct(v0, v0);
            double d01 = GetDotproduct(v0, v1);
            double d11 = GetDotproduct(v1, v1);
            double d20 = GetDotproduct(v2, v0);
            double d21 = GetDotproduct(v2, v1);
            double denom = d00 * d11 - d01 * d01;
            v = (d11 * d20 - d01 * d21) / denom;
            w = (d00 * d21 - d01 * d20) / denom;
            u = 1.0f - v - w;
        }
      
        static Vector3d Barycentric(const Vector3d& p, const Vector3d& a, const Vector3d& b, const Vector3d& c)
        {
            Vector3d v0 = GetMinus(b, a), v1 = GetMinus(c, a), v2 = GetMinus(p, a);
            double d00 = GetDotproduct(v0, v0);
            double d01 = GetDotproduct(v0, v1);
            double d11 = GetDotproduct(v1, v1);
            double d20 = GetDotproduct(v2, v0);
            double d21 = GetDotproduct(v2, v1);
            double denom = d00 * d11 - d01 * d01;
            double v = (d11 * d20 - d01 * d21) / denom;
            double w = (d00 * d21 - d01 * d20) / denom;
            double u = 1.0f - v - w;
            return Vector3d(u, v, w);
        }

       
        //Determine whether the three points are collinear
        template <class Type>
        static bool DetectColinear(const Type& v, const Type& s, const Type& e, const double& angle_match_error, const double& dis_match_error)
        {
            if (DetectCoincident(v, s, dis_match_error) || DetectCoincident(v, e, dis_match_error)) return true;//有其中两个点重合说明共线
            double angle = GetAngleBetween(v - s, e - s);//不共线就计算夹角，看是否为0或π
            if (IsAlmostZero_Double(angle, angle_match_error))return true;
            if (IsAlmostZero_Double(angle - Math_PI, angle_match_error))return true;
            return false;
        }
        //Determine whether two vectors are perpendicular
        template <class Type>
        static bool DetectVertical(const Type& direction_0, const Type& direction_1, const double& angle_match_error, const double& dis_match_error)
        {
            auto angle = GetAngleBetween(direction_0, direction_1);
            return IsAlmostZero_Double(angle - Math_PI / 2.0, angle_match_error);
        };

        template <class Type>
        static bool DetectVertical(const Type& seg_0_s, const Type& seg_0_e, const Type& seg_1_s, const Type& seg_1_e, const double& angle_match_error, const double& dis_match_error)
        {
            return DetectVertical(seg_0_e - seg_0_s, seg_1_e - seg_1_s, angle_match_error, dis_match_error);
        };

        template <class Type>
        static bool DetectVertical(const std::pair<Type, Type>& seg_0, const std::pair<Type, Type>& seg_1,
            const double& angle_match_error, const double& dis_match_error)
        {
            return DetectVertical(seg_0.first, seg_0.second, seg_1.first, seg_1.second, angle_match_error, dis_match_error);
        };
        //Determine whether two points coincide
        template <class Type>
        static bool DetectCoincident(const Type& v0, const Type& v1, const double& EPSILON = DOUBLE_EPSILON)
        {
            return IsAlmostZero_Double(GetDistance(v0, v1), EPSILON);
        }
        //Determine whether two planes are coplanar
        static bool DetectCoplanar(const Vector3d& planar_location_0, const Vector3d& planar_direction_0,
            const Vector3d& planar_location_1, const Vector3d& planar_direction_1,
            const double& angle_match_error, const double& dis_match_error)
        {
            auto angle = GetAngleBetween(planar_direction_0, planar_direction_1);
            if (IsAlmostZero_Double(angle - Math_PI, angle_match_error) || IsAlmostZero_Double(angle, angle_match_error))//Two parallel planes
            {
                double dis = GetLength(PlaneProject(planar_location_0, planar_direction_0, planar_location_1), planar_location_1);//Distance between parallel planes
                return IsAlmostZero_Double(dis, dis_match_error);
            }
            else
                return false;
        }
        //Determine whether two vectors are parallel
        template <class Type>
        static bool DetectParallel(const Type& direction_0, const Type& direction_1, const double& angle_match_error, const double& dis_match_error)
        {
            auto angle = GetAngleBetween(direction_0, direction_1);
            return (IsAlmostZero_Double(angle - Math_PI, angle_match_error) || IsAlmostZero_Double(angle, angle_match_error));
        };

        template <class Type>
        static bool DetectParallel(const std::pair<Type, Type>& seg_0, const std::pair<Type, Type>& seg_1,
            const double& angle_match_error, const double& dis_match_error)
        {
            return DetectParallel(seg_0.second - seg_0.first, seg_1.second - seg_1.first, angle_match_error, dis_match_error);
        };
        //Determine if the vectors are in the same direction
        template <class Type>
        static bool DetectCoDirection(const Type& direction_0, const Type& direction_1, const double& angle_match_error, const double& dis_match_error)
        {
            auto angle = GetAngleBetween(direction_0, direction_1);
            return (IsAlmostZero_Double(angle, angle_match_error));
        };
        // Determine whether two straight lines are collinear
        template <class Type>
        static bool DetectColinearDirection(const Type& location_0, const Type& direction_0, const Type& location_1, const Type& direction_1, const double& angle_match_error, const double& dis_match_error)
        {
            if (DetectParallel(direction_0, direction_1, angle_match_error, dis_match_error))//parallel
            {
                if (IsAlmostZero_Double(GetLength(location_0 - location_1), dis_match_error))//Starting point  coincident 
                    return true;
                return DetectParallel(direction_0, location_1 - location_0, angle_match_error, dis_match_error);//starting point not coincident but collinear
            }
            else
                return false;
        };
        //Determine whether two straight lines in two-dimensional space coincide
        //s_ 0, s_ 1 is the starting point of two straight lines, e_ 0, e_ 1 is the endpoint of two straight lines
        static bool DetectAlign2D(const Vector2d& s_0, const Vector2d& e_0, const Vector2d& s_1, const Vector2d& e_1, const double& angle_match_error, const double& dis_match_error)
        {
            double d0 = GetDistance(s_0, s_1);
            double d1 = GetDistance(s_0, e_1);
            double d2 = GetDistance(e_0, s_1);
            double d3 = GetDistance(e_0, e_1);
            if (IsAlmostZero_Double(d0, dis_match_error) && IsAlmostZero_Double(d3, dis_match_error))
                return true;
            if (IsAlmostZero_Double(d1, dis_match_error) && IsAlmostZero_Double(d2, dis_match_error))
                return true;
            return false;
        };
        //Determine whether two straight lines in three-dimensional space coincide
        static bool DetectAlign3D(const Vector3d& s_0, const Vector3d& e_0, const Vector3d& s_1, const Vector3d& e_1, const double& angle_match_error, const double& dis_match_error)
        {
            double d0 = GetDistance(s_0, s_1);
            double d1 = GetDistance(s_0, e_1);
            double d2 = GetDistance(e_0, s_1);
            double d3 = GetDistance(e_0, e_1);
            if (IsAlmostZero_Double(d0, dis_match_error) && IsAlmostZero_Double(d3, dis_match_error))
                return true;
            if (IsAlmostZero_Double(d1, dis_match_error) && IsAlmostZero_Double(d2, dis_match_error))
                return true;
            return false;
        };

        //this function has bug ;
        //Do not use it
        template <class Type>
        static bool DetectColinearSegment(const Type& s_0, const Type& e_0, const Type& s_1, const Type& e_1, const double& angle_match_error, const double& dis_match_error)
        {
            return DetectColinear_Direction(s_0, e_0 - s_0, s_1, e_1 - s_1, angle_match_error, dis_match_error);
        };

        //Clear duplicate and collinear points from the polygon, vecs is the point set of the polygon (in order)
        static std::vector<Vector2d> PolygonClear(const std::vector<Vector2d>& vecs,
            const double& angle_match_error, const double& dis_match_error)
        {
            //remove duplicate points
            std::vector<Vector2d> vecs_0;//Store point sets after deduplication
            for (int i = 0; i < vecs.size(); i++)
            {
                if (i != vecs.size() - 1)
                {
                    if (vecs_0.empty()) vecs_0.emplace_back(vecs[i]);
                    else
                    {
                        if (dis_match_error > 0)
                        {
                            if (!IsAlmostZero_Double(GetLength(vecs_0.back(), vecs[i]), dis_match_error))
                                vecs_0.emplace_back(vecs[i]);
                        }
                        else
                        {
                            if (!IsAlmostZero(GetLength(vecs_0.back(), vecs[i])))
                                vecs_0.emplace_back(vecs[i]);
                        }
                    }
                }
                else
                {
                    if (vecs_0.empty())
                        vecs_0.emplace_back(vecs[i]);
                    else
                    {
                        if (dis_match_error > 0)
                        {
                            if (!IsAlmostZero_Double(GetLength(vecs_0.back(), vecs[i]), dis_match_error) &&
                                !IsAlmostZero_Double(GetLength(vecs_0.front(), vecs[i]), dis_match_error))
                                vecs_0.emplace_back(vecs[i]);
                        }
                        else
                        {
                            if (!IsAlmostZero(GetLength(vecs_0.back(), vecs[i])) &&
                                !IsAlmostZero(GetLength(vecs_0.front(), vecs[i])))
                                vecs_0.emplace_back(vecs[i]);
                        }

                    }
                }
            }
            //remove collinear points
            std::vector<Vector2d> vecs_1;//Store point sets after removing duplicates and collinear points
            for (int i = 0; i < vecs_0.size(); i++)
            {
                auto pre_v = vecs_0[(i + vecs_0.size() - 1) % vecs_0.size()];//Previous vertex
                auto cur_v = vecs_0[i];//current vertex
                auto next_v = vecs_0[(static_cast<int64_t>(i) + 1) % vecs_0.size()];//Next vertex
                double angle = GetAngleBetween(cur_v - pre_v, next_v - cur_v);//An angle of 0 indicates that three points are collinear and cur can be removed_ V
                if (angle_match_error > 0.0)
                {
                    if (!IsAlmostZero_Double(angle, angle_match_error))
                        vecs_1.emplace_back(cur_v);
                }
                else
                {
                    if (!IsAlmostZero(angle))
                        vecs_1.emplace_back(cur_v);
                }
            }

            return vecs_1;
        };
        //Obtain Unit Square
        static std::vector<Vector2d> GetUnitSquare()
        {
            std::vector<Vector2d> square;
            square.push_back(Vector2d(-0.5, -0.5));
            square.push_back(Vector2d(0.5, -0.5));
            square.push_back(Vector2d(0.5, 0.5));
            square.push_back(Vector2d(-0.5, 0.5));
            return square;
        };
        //Obtaining the edge set of a cube
      
        static std::vector<std::pair<Vector3d, Vector3d>> GetUnitCubeFrame(const double& scale = 1.0)
        {
            Vector3d1 cube_vecs;//Unit cube

            cube_vecs.push_back(Vector3d(0.5, 0.5, 0.5));
            cube_vecs.push_back(Vector3d(-0.5, 0.5, 0.5));
            cube_vecs.push_back(Vector3d(-0.5, 0.5, -0.5));
            cube_vecs.push_back(Vector3d(0.5, 0.5, -0.5));

            cube_vecs.push_back(Vector3d(0.5, -0.5, 0.5));
            cube_vecs.push_back(Vector3d(-0.5, -0.5, 0.5));
            cube_vecs.push_back(Vector3d(-0.5, -0.5, -0.5));
            cube_vecs.push_back(Vector3d(0.5, -0.5, -0.5));

            auto sm = Functs::ScaleMatrix(Vector3d(scale, scale, scale));//Obtain transformation matrix
            cube_vecs = Functs::PosApplyMatrix(cube_vecs, sm);//coordinate transformation

            VectorPI1 frames_indexes;//Store the index for each edge
            frames_indexes.push_back(std::pair<int, int>(0, 1));
            frames_indexes.push_back(std::pair<int, int>(1, 2));
            frames_indexes.push_back(std::pair<int, int>(2, 3));
            frames_indexes.push_back(std::pair<int, int>(3, 0));
            frames_indexes.push_back(std::pair<int, int>(5, 4));
            frames_indexes.push_back(std::pair<int, int>(4, 7));
            frames_indexes.push_back(std::pair<int, int>(7, 6));
            frames_indexes.push_back(std::pair<int, int>(6, 5));
            frames_indexes.push_back(std::pair<int, int>(5, 1));
            frames_indexes.push_back(std::pair<int, int>(4, 0));
            frames_indexes.push_back(std::pair<int, int>(7, 3));
            frames_indexes.push_back(std::pair<int, int>(6, 2));

            std::vector<std::pair<Vector3d, Vector3d>> frames;//Set of edges
            for (auto& fi : frames_indexes)
                frames.push_back(std::pair<Vector3d, Vector3d>(cube_vecs[fi.first], cube_vecs[fi.second]));
            return frames;
        }
        //Divide each face of the scaled unit cube into triangular faces（get LiblgpTriMesh）
        //cube_vecs:the Point Set of a Cube
        static void GetUnitCube(Vector3d1& cube_vecs, Vector1i1& cube_face_id_0, Vector1i1& cube_face_id_1, Vector1i1& cube_face_id_2, const double& scale = 1.0)
        {
            cube_vecs.clear();
            cube_face_id_0.clear();
            cube_face_id_1.clear();
            cube_face_id_2.clear();

            cube_vecs.push_back(Vector3d(0.5, 0.5, 0.5));
            cube_vecs.push_back(Vector3d(-0.5, 0.5, 0.5));
            cube_vecs.push_back(Vector3d(-0.5, 0.5, -0.5));
            cube_vecs.push_back(Vector3d(0.5, 0.5, -0.5));

            cube_vecs.push_back(Vector3d(0.5, -0.5, 0.5));
            cube_vecs.push_back(Vector3d(-0.5, -0.5, 0.5));
            cube_vecs.push_back(Vector3d(-0.5, -0.5, -0.5));
            cube_vecs.push_back(Vector3d(0.5, -0.5, -0.5));

            auto sm = Functs::ScaleMatrix(Vector3d(scale, scale, scale));
            cube_vecs = Functs::PosApplyMatrix(cube_vecs, sm);

            Vector1i2 quad_faces;//Store vertex indices for each face
            quad_faces.push_back(Vector1i1{ 0, 1, 2, 3 });
            quad_faces.push_back(Vector1i1{ 5, 1, 0, 4 });
            quad_faces.push_back(Vector1i1{ 4, 0, 3, 7 });
            quad_faces.push_back(Vector1i1{ 5, 4, 7, 6 });
            quad_faces.push_back(Vector1i1{ 7, 3, 2, 6 });
            quad_faces.push_back(Vector1i1{ 6, 2, 1, 5 });

            for (auto qf : quad_faces)
            {
                cube_face_id_0.push_back(qf[2]);
                cube_face_id_1.push_back(qf[1]);
                cube_face_id_2.push_back(qf[0]);
                cube_face_id_0.push_back(qf[0]);
                cube_face_id_1.push_back(qf[3]);
                cube_face_id_2.push_back(qf[2]);
            }
        };

        static LiblgpTriMesh GetUnitCube(const double& scale = 1.0)
        {
            LiblgpTriMesh tm;
            GetUnitCube(tm.vecs, tm.face_id_0, tm.face_id_1, tm.face_id_2, scale);
            return tm;
        }
        //Returns the direction vector representing the direction of the coordinate axis
        static Vector3d1 GetAxisAlignDirections()
        {
            Vector3d1 directions;
            directions.push_back(Vector3d(-1, 0, 0));
            directions.push_back(Vector3d(1, 0, 0));
            directions.push_back(Vector3d(0, -1, 0));
            directions.push_back(Vector3d(0, 1, 0));
            directions.push_back(Vector3d(0, 0, -1));
            directions.push_back(Vector3d(0, 0, 1));
            return directions;
        }
        //Obtain a list of vectors containing 24 rotated Euler angles
        static Vector3d1 EmumerateRotations()
        {
            Vector3d1 rotations;
            rotations.emplace_back(Vector3d(90.0, 90.0, 180.0));
            rotations.emplace_back(Vector3d(-90.0, 0.0, 90.0));
            rotations.emplace_back(Vector3d(0.0, 180.0, 0.0));

            rotations.emplace_back(Vector3d(90.0, 0.0, 180.0));
            rotations.emplace_back(Vector3d(0.0, 0.0, 90.0));
            rotations.emplace_back(Vector3d(0.0, -90.0, 0.0));

            rotations.emplace_back(Vector3d(0.0, 0.0, 0.0));
            rotations.emplace_back(Vector3d(90.0, 0.0, 90.0));
            rotations.emplace_back(Vector3d(90.0, -90.0, 180.0));

            rotations.emplace_back(Vector3d(0.0, 90.0, 0.0));
            rotations.emplace_back(Vector3d(180.0, 0.0, 90.0));
            rotations.emplace_back(Vector3d(90.0, -180.0, 180.0));

            rotations.emplace_back(Vector3d(0.0, 180.0, 90.0));
            rotations.emplace_back(Vector3d(-90.0, 90.0, 90.0));
            rotations.emplace_back(Vector3d(90.0, 180.0, 0.0));

            rotations.emplace_back(Vector3d(0.0, 0.0, 180.0));
            rotations.emplace_back(Vector3d(0.0, -90.0, 90.0));
            rotations.emplace_back(Vector3d(90.0, 0.0, -90.0));

            rotations.emplace_back(Vector3d(-90.0, 0.0, -90.0));
            rotations.emplace_back(Vector3d(180.0, 90.0, 90.0));
            rotations.emplace_back(Vector3d(0.0, -180.0, 180.0));

            rotations.emplace_back(Vector3d(90.0, 0.0, 0.0));
            rotations.emplace_back(Vector3d(90.0, -90.0, 90.0));
            rotations.emplace_back(Vector3d(0.0, 0.0, 270.0));

            for (auto& rotation : rotations)//Convert angles to radians
            {
                rotation[0] = rotation[0] / 180.0 * Math_PI;
                rotation[1] = rotation[1] / 180.0 * Math_PI;
                rotation[2] = rotation[2] / 180.0 * Math_PI;
            }
            return rotations;
        };


#pragma endregion

#pragma region BasicMathFunctions
        //Returns a random double precision floating-point number within a specified range
        static double RandomDouble(const double& min_d = 0.0, const double& max_d = 1.0)
        {
            std::uniform_real_distribution<> dis(min_d, max_d);//Uniformly distributed random numbers，range：[min_d,max_d]
            return dis(MATHGEN);
        }

        //select a number among 0,1,2,...,s-1
        //s should be positive int
        static int RandomInt(const int& s)
        {
            if (s == 1) return 0;
            double x = 1.0 / s;
            double d = RandomDouble();

            for (int i = 0; i < s; i++)
            {
                double min_d = i * x;
                double max_d = (static_cast<int64_t>(i) + 1)* x;
                if (d >= min_d && d <= max_d)
                {
                    return i;
                }
            }

            return -1;
        }

        /*static double RandomDD(const double& min_d = 0.0, const double& max_d = 1.0)//redundancy？
        {
            std::uniform_real_distribution<> dis(min_d, max_d);
            return dis(MATHGEN);
        }*/

        //select a number among 0,1,2,...,s-1
        //s should be positive int
        /*static int RandomII(int s)
        {
            if (s == 1) return 0;
            double x = 1.0 / s;
            double d = RandomDD();

            for (int i = 0; i < s; i++)
            {
                double min_d = i * x;
                double max_d = (static_cast<int64_t>(i) + 1)* x;
                if (d >= min_d && d <= max_d)
                {
                    return i;
                }
            }

            return -1;
        }*/

        //Calculate vector cross product
        static Vector3d GetCrossproduct(const Vector3d& v1, const Vector3d& v2) {
            return glm::cross(v1, v2);
        }
        //Calculate vector dot product
        template <class Type>
        static double GetDotproduct(const Type& v1, const Type& v2) {
            return glm::dot(v1, v2);
        }
        //Calculate difference
        template <class Type>
        static Type GetMinus(const Type& a, const Type& b)
        {
            Type c = a;
            for (int i = 0; i < a.length(); i++)
                c[i] = a[i] - b[i];
            return c;
        }
        //Judging equality
        static bool AreAlmostEqual(const double& value1, const double& value2) {
            if (value1 == value2) {
                return true;
            }
            double eps = (glm::abs(value1) + glm::abs(value2) + 10.0) * DOUBLE_EPSILON;
            double delta = value1 - value2;
            return (-eps < delta) && (eps > delta);
        }

        static bool AreAlmostEqual_Double(const double& value1, const double& value2, const double& EPSILON) {
            return IsAlmostZero_Double(value1 - value2, EPSILON);
        }
        //Returns true if the number of a double type is approximately zero
        static bool IsAlmostZero(const double& value) {
            return (value < DOUBLE_EPSILON) && (value > -DOUBLE_EPSILON);
        }
        static bool IsAlmostZero_Double(const double& value, const double& EPSILON) {
            return (value < EPSILON) && (value > -EPSILON);
        }

        /// Returns true if two given floating point numbers are epsilon-equal.
        /// Method automatically adjust the epsilon to the absolute size of given numbers.
        static bool AreAlmostEqual(const float& value1, const float& value2) {
            // in case they are Infinities (then epsilon check does not work)
            if (value1 == value2) {
                return true;
            }
            // computes (|value1-value2| / (|value1| + |value2| + 10.0)) < SINGLE_EPSILON
            float eps = (float)((glm::abs(value1) + glm::abs(value2) + 10.0) * SINGLE_EPSILON);
            float delta = value1 - value2;
            return (-eps < delta) && (eps > delta);
        }
        //Set the parts of a three-dimensional space point that are approximately zero in the three dimensions to zero
        static void ZeroVector(Vector3d& v)
        {
            if (IsAlmostZero(v[0]))v[0] = 0.0;
            if (IsAlmostZero(v[1]))v[1] = 0.0;
            if (IsAlmostZero(v[2]))v[2] = 0.0;
        }
        //Calculate maximum value
        template<class Type>
        static double GetMax(const Type& vec)
        {
            double maxd = vec[0];
            for (int i = 0; i < vec.length(); i++)
                maxd = max(maxd, vec[i]);
            return maxd;
        }
        //Calculate the maximum and minimum values of a one-dimensional array
        static void GetMinMax(const Vector1d1& ds, double& mind, double& maxd)
        {
            mind = ds[0];
            maxd = ds[0];
            for (auto& d : ds)
            {
                mind = min(mind, d);
                maxd = max(maxd, d);
            }
        }
        //Calculate the maximum and minimum values of a two-dimensional array
        static void GetMinMax(const Vector1d2& ds, double& mind, double& maxd)
        {
            mind = ds.front().front();
            maxd = ds.front().front();
            for (auto& d : ds)
            {
                for (auto& d_ : d)
                {
                    mind = min(mind, d_);
                    maxd = max(maxd, d_);
                }
            }
        }
        //Insert a non repeating element into the array, and return false if the element already exists in the array.Otherwise, insert and return true
        template<class Type>
        static bool VectorInsertNoDuplicate(std::vector<Type>& vecs, const Type& element)
        {
            if (CheckContain(vecs, element))
                return false;

            vecs.emplace_back(element);
            return true;
        }
       
        template<class Type>
        static void VectorInsertNoDuplicate(std::vector<Type>& vecs, const std::vector<Type>& elements)
        {
            for (auto& element : elements)
                VectorInsertNoDuplicate(vecs, element);
        }
        //Merge two arrays
        template<class Type>
        static std::vector<Type> VectorMerge(const std::vector<Type>& vecs_0, const std::vector<Type>& vecs_1)
        {
            std::vector<Type> result = vecs_0;
            result.insert(result.end(), vecs_1.begin(), vecs_1.end());
            return result;
        }
        //Check if the array contains the given element
        template<class Type>
        static bool CheckContain(const Vector1<Type>& vecs, const Type& element)
        {
            return std::find(vecs.begin(), vecs.end(), element) != vecs.end();
        }

        template<class Type>
        static bool CheckContain(const Vector2<Type>& vecs, const Type& element)
        {
            return VectorIndex(vecs, element) >= 0;
            //return std::find(vecs.begin(), vecs.end(), element) != vecs.end();
        }

        //Return the index of the element in an array,or -1 if there is no such element
        template <class Type>
        static int VectorIndex(const std::vector <Type>& vecs, const Type& element)
        {
            for (int i = 0; i < vecs.size(); i++)
            {
                if (vecs[i] == element)
                {
                    return i;
                }
            }
            return -1;
        }

        template <class Type>
        static int VectorIndex(const std::vector <std::vector <Type>>& vecs, const Type& element)
        {
            for (int i = 0; i < vecs.size(); i++)
            {
                if (VectorIndex(vecs[i], element) >= 0)
                    return i;
            }
            return -1;
        }

        //Return the size of a two-dimensional array
        template <class Type>
        static int VectorSize(const Vector2<Type>& vecs)
        {
            int nb = 0;
            for (auto& vec : vecs)
                nb += vec.size();
            return nb;
        }
        //Return the size of a three-dimensional array
        template <class Type>
        static int VectorSize(const Vector3<Type>& vecs)
        {
            int nb = 0;
            for (auto& vec : vecs)
                nb += VectorSize(vec);
            return nb;
        }

      
        //Add a given value to each element in the array
        template <class Type>
        static std::vector <Type> VectorAdd(const std::vector <Type>& vecs, const Type& element)
        {
            std::vector <Type> result = vecs;
            for (auto& r : result)
                r += element;
            return result;
        }
        //Return true if the array contains the given element
        static bool VectorContain(const std::vector<int>& vecs, const int& element)
        {
            for (int i = 0; i < vecs.size(); i++)
            {
                if (vecs[i] == element)
                    return true;
            }

            return false;
        }
        //Return the index of the given element in the array, or -1 if it does not contain the element
        static int VectorContainReturnIndex(const std::vector<int>& vecs, const int& element)
        {
            for (int i = 0; i < vecs.size(); i++)
            {
                if (vecs[i] == element)
                    return i;
            }

            return -1;
        }
        //Return true if an array in two-dimensional space contains a given element (ordered)
        static bool VectorContainForSpecialCase(const std::vector<std::vector<int>>& vecs, const std::vector<int>& element)
        {
            for (int i = 0; i < vecs.size(); i++)
            {
                if (vecs[i][0] == element[0] && vecs[i][1] == element[1])
                    return true;
            }

            return false;
        }
        //Return true if an array in two-dimensional space contains a given element (unordered)
        static bool VectorContainForSpecialCase1(const std::vector<std::vector<int>>& vecs, const std::vector<int>& element)
        {
            for (int i = 0; i < vecs.size(); i++)
            {
                if (vecs[i][0] == element[0] && vecs[i][1] == element[1]) return true;
                if (vecs[i][0] == element[1] && vecs[i][1] == element[0]) return true;
            }

            return false;
        }
        //Return the index of the given element in an array in two-dimensional space, or -1 if there is no such element(unordered)
        static int VectorContainForSpecialCase2(const std::vector<std::vector<int>>& vecs, const int& element_0, const int& element_1)
        {
            for (int i = 0; i < vecs.size(); i++)
            {
                if (vecs[i][0] == element_0 && vecs[i][1] == element_1) return i;
                if (vecs[i][0] == element_1 && vecs[i][1] == element_0) return i;
            }
            return -1;
        }
        //Return the index of the given element in an array in two-dimensional space, or -1 if there is no such element(ordered)
        static int VectorContainForSpecialCase3(const std::vector<std::vector<int>>& vecs, const int& element_0, const int& element_1)
        {
            for (int i = 0; i < vecs.size(); i++)
            {
                if (vecs[i][0] == element_0 && vecs[i][1] == element_1) return i;
            }
            return -1;
        }
#pragma endregion


#pragma region Transformation
        //Reduce 3d to 2d, remove the z-axis coordinates
        static Vector2d Vector3dto2d(const Vector3d& v)
        {
            return Vector2d(v[0], v[1]);
        }
        //Reduce  3d array to  2d array (array dimension is one-dimensional)
        static Vector2d1 Vector3dto2d(const Vector3d1& vecs_3d)
        {
            Vector2d1 vecs_2d;
            for (auto& v : vecs_3d)
                vecs_2d.emplace_back(Vector3dto2d(v));
            return vecs_2d;
        }
        //Reduce  3d array to  2d array (array dimension is two-dimensional)
        static Vector2d2 Vector3dto2d(const Vector3d2& vecs_3d)
        {
            Vector2d2 vecs_2d;
            for (auto v : vecs_3d)
                vecs_2d.emplace_back(Vector3dto2d(v));
            return vecs_2d;
        }
        //2d to 3d. z:set the value of the z-axis
        static Vector3d Vector2dto3d(const Vector2d& v, const double& z = 0.0)
        {
            return Vector3d(v[0], v[1], z);
        }
        //2d array to 3d array (array dimension is one-dimensional)
        static Vector3d1 Vector2dto3d(const Vector2d1& vecs_2d, double z = 0.0)
        {
            Vector3d1 vecs_3d;
            for (auto v : vecs_2d)
                vecs_3d.emplace_back(Vector2dto3d(v, z));
            return vecs_3d;
        }
        //2d array to 3d array (array dimension is two-dimensional)
        static Vector3d2 Vector2dto3d(const Vector2d2& vecs_2d, double z = 0.0)
        {
            Vector3d2 vecs_3d;
            for (auto v : vecs_2d)
                vecs_3d.emplace_back(Vector2dto3d(v, z));
            return vecs_3d;
        }
        //Apply a 3d vector v to a 4x4 transformation matrix M and return the transformed result
        static Vector3d VecApplyMatrix(const Vector3d& v, const  glm::dmat4& M)
        {
            return Vector3d(M * glm::vec4(v, 1.0)) - Vector3d(M * glm::vec4(Vector3d(0.0, 0.0, 0.0), 1.0));
        }
        //Affine transformation of one-dimensional vector array
        static Vector3d1 VecApplyMatrix(const Vector3d1& vecs, const glm::dmat4& M)
        {
            Vector3d1 ps;
            for (auto& p : vecs)
                ps.emplace_back(VecApplyMatrix(p, M));
            return ps;
        }
        //Affine transformation of  two-dimensional vector array
        static Vector3d2 VecApplyMatrix(const Vector3d2& veces, const glm::dmat4& M)
        {
            Vector3d2 pses;
            for (auto vecs : veces)
                pses.emplace_back(VecApplyMatrix(vecs, M));
            return pses;
        }
        // Coordinate transformation of point
        static Vector3d PosApplyMatrix(const Vector3d& v, const glm::dmat4& M)
        {
            return Vector3d(M * glm::vec4(v, 1.0));
        }
        //Coordinate transformation of one-dimensional point array
        static Vector3d1 PosApplyMatrix(const Vector3d1& vecs, const glm::dmat4& M)
        {
            Vector3d1 ps;
            for (auto& p : vecs)
                ps.emplace_back(PosApplyMatrix(p, M));
            return ps;
        }

        
        //Coordinate transformation of two-dimensional point array
        static Vector3d2 PosApplyMatrix(const Vector3d2& veces, const glm::dmat4& M)
        {
            Vector3d2 pses;
            for (auto vecs : veces)
                pses.emplace_back(PosApplyMatrix(vecs, M));
            return pses;
        }
        //Coordinate transformation of three-dimensional point array
        static Vector3d3 PosApplyMatrix(const Vector3d3& veces, const glm::dmat4& M)
        {
            Vector3d3 pses;
            for (auto vecs : veces)
                pses.emplace_back(PosApplyMatrix(vecs, M));
            return pses;
        }
        // Coordinate transformation of point.form like (point,point)
        static std::pair<Vector3d, Vector3d> PosApplyMatrix(const std::pair<Vector3d, Vector3d>& vecs, const glm::dmat4& M)
        {
            return std::pair<Vector3d, Vector3d>(PosApplyMatrix(vecs.first, M), PosApplyMatrix(vecs.second, M));
        }
        //Get the rotation matrix
        //x[0] y[0] z[0] 0
        //x[1] y[1] z[1] 0
        //x[2] y[2] z[2] 0
        // 0    0    0   1
        static glm::dmat4 RotationMatrixXYZ(const Vector3d& xx, const Vector3d& yy, const Vector3d& zz)
        {
            auto x = xx;
            auto y = yy;
            auto z = zz;

            ZeroVector(x);
            ZeroVector(y);
            ZeroVector(z);
            x = x / (double)GetLength(x);
            y = y / (double)GetLength(y);
            z = z / (double)GetLength(z);

            glm::dmat4  rotationMatrix;

            rotationMatrix[0][0] = x[0];
            rotationMatrix[0][1] = y[0];
            rotationMatrix[0][2] = z[0];
            rotationMatrix[0][3] = 0.0;

            rotationMatrix[1][0] = x[1];
            rotationMatrix[1][1] = y[1];
            rotationMatrix[1][2] = z[1];
            rotationMatrix[1][3] = 0.0;

            rotationMatrix[2][0] = x[2];
            rotationMatrix[2][1] = y[2];
            rotationMatrix[2][2] = z[2];
            rotationMatrix[2][3] = 0.0;

            rotationMatrix[3][0] = 0.0;
            rotationMatrix[3][1] = 0.0;
            rotationMatrix[3][2] = 0.0;
            rotationMatrix[3][3] = 1.0;

            return rotationMatrix;

        }

        //Generate a rotation matrix, which will rotate vector o around the rotation axis n to the same direction as vector t
        static glm::dmat4 RotationMatrix(const Vector3d& o, const Vector3d& t, const Vector3d& n)
        {
            double angle = GetAngleBetween(o, t);

            if (IsAlmostZero(angle))
            {
                glm::dmat4  rotationMatrix;

                rotationMatrix[0][0] = 1.0;
                rotationMatrix[0][1] = 0.0;
                rotationMatrix[0][2] = 0.0;
                rotationMatrix[0][3] = 0.0;
                rotationMatrix[1][0] = 0.0;
                rotationMatrix[1][1] = 1.0;
                rotationMatrix[1][2] = 0.0;
                rotationMatrix[1][3] = 0.0;
                rotationMatrix[2][0] = 0.0;
                rotationMatrix[2][1] = 0.0;
                rotationMatrix[2][2] = 1.0;
                rotationMatrix[2][3] = 0.0;
                rotationMatrix[3][0] = 0.0;
                rotationMatrix[3][1] = 0.0;
                rotationMatrix[3][2] = 0.0;
                rotationMatrix[3][3] = 1.0;

                return rotationMatrix;
            }
            else
            {
                return RotationMatrix(n, angle);
            }
        }
        /*//Adjust the value of the rotation axis n based on the given vector v
        static void AAA(const Vector3d& v, Vector3d& n)
        {
            auto a = v[0];
            auto b = v[1];
            auto c = v[2];
            bool bx = IsAlmostZero(a);
            bool by = IsAlmostZero(b);
            bool bz = IsAlmostZero(c);

            if (bx && by && bz)
            {
                std::cerr << "if (bx&&by&&bz)" << std::endl;
                system("pause");
            }

            if (bx && by && !bz)
            {
                n[0] = 1.0;
                n[1] = 1.0;
                n[2] = 0.0;
            }
            if (bx && !by && bz)
            {
                n[0] = 1.0;
                n[1] = 0.0;
                n[2] = 1.0;
            }

            if (bx && !by && !bz)
            {
                n[0] = 1.0;
                n[1] = 1.0;
                n[2] = -b / c;
            }

            if (!bx && by && bz)
            {
                n[0] = 0.0;
                n[1] = 1.0;
                n[2] = 1.0;
            }

            if (!bx && by && !bz)
            {
                n[0] = 1.0;
                n[1] = 1.0;
                n[2] = -a / c;
            }
            if (!bx && !by && bz)
            {
                n[0] = 1.0;
                n[1] = -a / b;
                n[2] = 1.0;
            }

            if (!bx && !by && !bz)
            {
                n[0] = 1.0;
                n[1] = 1.0;
                n[2] = -(a + b) / c;
            }

        }*/
        //Generate a rotation matrix based on the starting vector o and the target vector t, which is used to rotate vector o to the same direction as vector t
        static glm::dmat4 RotationMatrix(const Vector3d& o, const Vector3d& t)
        {
            Vector3d n = GetCrossproduct(o, t);
            double angle = GetAngleBetween(o, t);

           /* if (IsAlmostZero(angle - Math_PI))
            {
                AAA(o, n);
            }*/

            if (IsAlmostZero(angle))
            {
                glm::dmat4  rotationMatrix;

                rotationMatrix[0][0] = 1.0;
                rotationMatrix[0][1] = 0.0;
                rotationMatrix[0][2] = 0.0;
                rotationMatrix[0][3] = 0.0;
                rotationMatrix[1][0] = 0.0;
                rotationMatrix[1][1] = 1.0;
                rotationMatrix[1][2] = 0.0;
                rotationMatrix[1][3] = 0.0;
                rotationMatrix[2][0] = 0.0;
                rotationMatrix[2][1] = 0.0;
                rotationMatrix[2][2] = 1.0;
                rotationMatrix[2][3] = 0.0;
                rotationMatrix[3][0] = 0.0;
                rotationMatrix[3][1] = 0.0;
                rotationMatrix[3][2] = 0.0;
                rotationMatrix[3][3] = 1.0;

                return rotationMatrix;
            }
            else
            {
                return RotationMatrix(n, angle);
            }
        }


        //Generate a rotation matrix with a given axis n as the rotation axis and an angle of rotation
        static glm::dmat4 RotationMatrix(const Vector3d& n, const double& angle)
        {
            //return glm::rotate(angle, n);
            double u = n[0];
            double v = n[1];
            double w = n[2];

            glm::dmat4  rotationMatrix;

            double L = (u * u + v * v + w * w);

            //angle = angle * M_PI / 180.0; //converting to radian value
            double u2 = u * u;
            double v2 = v * v;
            double w2 = w * w;

            rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
            rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
            rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
            rotationMatrix[0][3] = 0.0;

            rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
            rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
            rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
            rotationMatrix[1][3] = 0.0;

            rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
            rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
            rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
            rotationMatrix[2][3] = 0.0;

            rotationMatrix[3][0] = 0.0;
            rotationMatrix[3][1] = 0.0;
            rotationMatrix[3][2] = 0.0;
            rotationMatrix[3][3] = 1.0;

            return rotationMatrix;
        }
       
        //Generate a translation matrix
        static glm::dmat4 TranslationMatrix(const Vector3d& v)
        {
            glm::dmat4  translationMatrix;
            translationMatrix[0][0] = 1.0;
            translationMatrix[0][1] = 0.0;
            translationMatrix[0][2] = 0.0;
            translationMatrix[0][3] = 0.0;

            translationMatrix[1][0] = 0.0;
            translationMatrix[1][1] = 1.0;
            translationMatrix[1][2] = 0.0;
            translationMatrix[1][3] = 0.0;

            translationMatrix[2][0] = 0.0;
            translationMatrix[2][1] = 0.0;
            translationMatrix[2][2] = 1.0;
            translationMatrix[2][3] = 0.0;

            translationMatrix[3][0] = v[0];
            translationMatrix[3][1] = v[1];
            translationMatrix[3][2] = v[2];
            translationMatrix[3][3] = 1.0;

            return translationMatrix;
        }
        //Generate a scaling matrix,  v is the scaling factor
        //v[0]  0   0   0
        //0    v[1] 0   0
        //0     0   v[2] 0
        //0    0    0    1
        static glm::dmat4 ScaleMatrix(const Vector3d& v)
        {
            glm::dmat4  translationMatrix;
            translationMatrix[0][0] = v[0];
            translationMatrix[0][1] = 0.0;
            translationMatrix[0][2] = 0.0;
            translationMatrix[0][3] = 0.0;

            translationMatrix[1][0] = 0.0;
            translationMatrix[1][1] = v[1];
            translationMatrix[1][2] = 0.0;
            translationMatrix[1][3] = 0.0;

            translationMatrix[2][0] = 0.0;
            translationMatrix[2][1] = 0.0;
            translationMatrix[2][2] = v[2];
            translationMatrix[2][3] = 0.0;

            translationMatrix[3][0] = 0;
            translationMatrix[3][1] = 0;
            translationMatrix[3][2] = 0;
            translationMatrix[3][3] = 1.0;

            return translationMatrix;
        }

        //Rotate around a fixed axis. p:the point to rotate, angle:the angle of rotation, n:the axis of rotation
        static Vector3d RotationAxis(const Vector3d& p, const double& angle, const Vector3d& n)
        {
            //auto m = RotationMatrix(n, angle);
            //return PosApplyM(p, m);

            auto rtv = glm::rotate(angle, n) * glm::dvec4(p, 1.0);
            return Vector3d(rtv[0], rtv[1], rtv[2]);

            /*
            glm::dmat4 inputMatrix(0.0);
            inputMatrix[0][0] = p[0];
            inputMatrix[1][0] = p[1];
            inputMatrix[2][0] = p[2];
            inputMatrix[3][0] = 1.0;
            double u = n[0];
            double v = n[1];
            double w = n[2];

            glm::dmat4  rotationMatrix;

            double L = (u * u + v * v + w * w);

            //angle = angle * M_PI / 180.0; //converting to radian value
            double u2 = u * u;
            double v2 = v * v;
            double w2 = w * w;

            rotationMatrix[0][0] = (u2 + (v2 + w2) * glm::cos(angle)) / L;
            rotationMatrix[0][1] = (u * v * (1 - glm::cos(angle)) - w * glm::sqrt(L) * glm::sin(angle)) / L;
            rotationMatrix[0][2] = (u * w * (1 - glm::cos(angle)) + v * glm::sqrt(L) * glm::sin(angle)) / L;
            rotationMatrix[0][3] = 0.0;

            rotationMatrix[1][0] = (u * v * (1 - glm::cos(angle)) + w * glm::sqrt(L) * glm::sin(angle)) / L;
            rotationMatrix[1][1] = (v2 + (u2 + w2) * glm::cos(angle)) / L;
            rotationMatrix[1][2] = (v * w * (1 - glm::cos(angle)) - u * glm::sqrt(L) * glm::sin(angle)) / L;
            rotationMatrix[1][3] = 0.0;

            rotationMatrix[2][0] = (u * w * (1 - glm::cos(angle)) - v * glm::sqrt(L) * glm::sin(angle)) / L;
            rotationMatrix[2][1] = (v * w * (1 - glm::cos(angle)) + u * glm::sqrt(L) * glm::sin(angle)) / L;
            rotationMatrix[2][2] = (w2 + (u2 + v2) * glm::cos(angle)) / L;
            rotationMatrix[2][3] = 0.0;

            rotationMatrix[3][0] = 0.0;
            rotationMatrix[3][1] = 0.0;
            rotationMatrix[3][2] = 0.0;
            rotationMatrix[3][3] = 1.0;

            double outputMatrix[4][1];

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 1; j++)
                {
                    outputMatrix[i][j] = 0;
                    for (int k = 0; k < 4; k++)
                    {
                        outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
                    }
                }
            }
            return Vector3d(outputMatrix[0][0], outputMatrix[1][0], outputMatrix[2][0]);
            */
        }
        //The 2d point p rotates around the specified center point and returns the rotated 2d point
        static Vector2d RotationAxis2d(const Vector2d& p, const double& angle, const Vector2d& center)
        {
            Vector3d r = RotationAxis(Vector3d(p[0] - center[0], 0.0, p[1] - center[1]),
                angle, Vector3d(0.0, 1.0, 0.0)) + Vector3d(center[0], 0.0, center[1]);
            return Vector2d(r[0], r[2]);
        }
        //The 3d point p rotates around the specified ray axis  and returns the rotated 3d point
        static Vector3d RotationAxis(const Vector3d& p, const double& angle, const Vector3d& ray_point, const Vector3d& ray_vector)
        {
            return RotationAxis(p - ray_point, angle, ray_vector) + ray_point;
        }
        //The 3d point array p rotates around the specified ray axis  and returns the rotated 3d point(array dimension:one-dimensional)
        static void RotationAxis(Vector3d1& points, const double& angle, const Vector3d& ray_point, const Vector3d& ray_vector)
        {
            for (int i = 0; i < points.size(); i++)
                points[i] = RotationAxis(points[i], angle, ray_point, ray_vector);
        }
        //The 3d point array p rotates around the specified ray axis  and returns the rotated 3d point(array dimension:two-dimensional)
        static void RotationAxis(Vector3d2& points, const double& angle, const Vector3d& ray_point, const Vector3d& ray_vector)
        {
            for (int i = 0; i < points.size(); i++)
                RotationAxis(points[i], angle, ray_point, ray_vector);
        }
        //The 3d point array p rotates around the specified ray axis  and returns the rotated 3d point(array dimension:three-dimensional)
        static void RotationAxis(Vector3d3& pointses, const double& angle, const Vector3d& ray_point, const Vector3d& ray_vector)
        {
            for (int i = 0; i < pointses.size(); i++)
                RotationAxis(pointses[i], angle, ray_point, ray_vector);
        }
        //Perform linear transformations on vectors
        static Vector3d Translate(const Vector3d& p, const Vector3d& v)
        {
            return p + v;
        }
        //Perform  linear transformation on  vector array(array dimension:one-dimensional)
        static void Translate(Vector3d1& points, const Vector3d& v)
        {
            for (int i = 0; i < points.size(); i++)
                points[i] = Translate(points[i], v);
        }
        //Perform  linear transformation on  vector array(array dimension:two-dimensional)
        static void Translate(Vector3d2& points, const Vector3d& v)
        {
            for (int i = 0; i < points.size(); i++)
                Translate(points[i], v);
        }
#pragma endregion


#pragma region IOFunctions
        //Check if the file with the specified path exists
        static bool LoadExisting(const std::string& path)
        {
            std::ifstream file(path, std::ios::in);
            if (!file) return false;
            return true;
        }
        //Load vector array data from a file at the given path（array dimension:three-dimensional )
        static bool LoadVectors(const std::string& path, Vector3d3& vec_3)
        {
            //zigzag_final_path
            int nb_0, nb_1, nb_2;//The number of 3d vectors, the number of 2d vectors in each 3d vector, the number of 1d vectors in each 2d vector
            std::ifstream file(path, std::ios::in);

            if (!file) return false;

            file >> nb_0;
            for (int i = 0; i < nb_0; i++)
            {
                file >> nb_1;
                Vector3d2 vec_2;
                for (int j = 0; j < nb_1; j++)
                {
                    file >> nb_2;
                    Vector3d1 vec_1(nb_2, Vector3d(0.0, 0.0, 0.0));
                    for (int k = 0; k < nb_2; k++)
                        file >> vec_1[k][0] >> vec_1[k][1] >> vec_1[k][2];
                    vec_2.emplace_back(vec_1);
                }
                vec_3.emplace_back(vec_2);
            }
            file.clear();
            file.close();

            return true;
        }
        //Load vector array data from a file at the given path（array dimension:one-dimensional )
        static bool LoadVectors(const std::string& path, Vector3d1& vec_3)
        {
            //zigzag_final_path
            std::ifstream file(path, std::ios::in);

            if (!file) return false;

            int nb;//The number of vectors
            file >> nb;
            for (int i = 0; i < nb; i++)
            {
                vec_3.emplace_back(Vector3d());
                file >> vec_3.back()[0] >> vec_3.back()[1] >> vec_3.back()[2];
            }
            file.clear();
            file.close();

            return true;
        }
        //Output the data of a vector array to a file at the given path（array dimension:three-dimensional )
        static void OutputVectors(const std::string& out_path, const Vector3d3& vecs)
        {
            std::ofstream file(out_path);
            file << vecs.size() << std::endl;

            for (int i = 0; i < vecs.size(); i++) {
                file << vecs[i].size() << std::endl;
                for (int j = 0; j < vecs[i].size(); j++) {
                    file << vecs[i][j].size() << std::endl;
                    for (int k = 0; k < vecs[i][j].size(); k++)
                        file << vecs[i][j][k][0] << " " << vecs[i][j][k][1] << " " << vecs[i][j][k][2] << " ";
                    file << "" << std::endl;
                }
            }

            file.clear();
            file.close();
        }
        //Output the data of a vector array to a file at the given path（array dimension:one-dimensional )
        static void OutputVectors(const std::string& out_path, const Vector3d1& vecs)
        {
            std::ofstream file(out_path);
            file << vecs.size() << std::endl;
            for (int i = 0; i < vecs.size(); i++)
                file << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
            file.clear();
            file.close();
        }

#if !defined(UNICODE) && !defined(_UNICODE) && !defined(__APPLE__)
        //Load a DLL file with the specified path
        static HMODULE LoadHMODULE(const string& dll_path)
        {
 
            if (!DetectExisting(dll_path))
            {
                std::string root_path = Functs::WinGetCurDirectory();

                std::string str;
                str += "The dll does not exist: " + dll_path + ";\n";
                str += "The current running directory: " + root_path + ";\n";
                str += "Please gurrentee the dll is in the right place;\n";
                MAssert(str);
            }

            HMODULE hModule = LoadLibrary(_T(dll_path.c_str()));
            if (!hModule)
            {
                DWORD dw = GetLastError(); // returns 0xc1 (193)
                MAssert("LoadLibrary failed with error code " + std::to_string(dw));
            }
            else
                std::cerr << "LoadLibrary success\n";

            return hModule;
        };

#endif
        //Load a 3D model file in. obj format
        //coords:store the vertex coordinates of the model(Every three consecutive elements represent a vertex coordinate)
        //tris:Store vertex index of triangular faces
        static void LoadObj3d(const char* path, std::vector<double>& coords, std::vector<int>& tris)
        {
            //Extract the first integer value from a string
            auto get_first_integer = [](const char* v)
            {
                int ival;
                std::string s(v);
                std::replace(s.begin(), s.end(), '/', ' ');
                sscanf(s.c_str(), "%d", &ival);//Read the first integers into ival
                return ival;
            };

            double x, y, z;
            char line[1024], v0[1024], v1[1024], v2[1024];

            // open the file, return if open fails
            FILE* fp = fopen(path, "r");
            if (!Functs::DetectExisting(path))
            {
                Functs::MAssert("This file does not exist: " + std::string(path));
                return;
            };

            while (fgets(line, 1024, fp))
            {
                if (line[0] == 'v')
                {
                    sscanf(line, "%*s%lf%lf%lf", &x, &y, &z);
                    coords.push_back(x);
                    coords.push_back(y);
                    coords.push_back(z);
                }
                else
                {
                    if (line[0] == 'f')
                    {
                        sscanf(line, "%*s%s%s%s", v0, v1, v2);
                        tris.push_back(get_first_integer(v0) - 1);
                        tris.push_back(get_first_integer(v1) - 1);
                        tris.push_back(get_first_integer(v2) - 1);
                    }
                }
            }
            fclose(fp);
        };
        //Load a 3D model file in. obj format
        //vecs:store vertex coordinates
        //face_id_0(1)(2):store the index of triangular faces
        static void LoadObj3d(const char* path_, Vector3d1& vecs, Vector1i1& face_id_0, Vector1i1& face_id_1, Vector1i1& face_id_2)
        {
            std::string path = path_;
            if (path.substr(path.size() - 3, path.size()) == "obj")
            {
                std::vector<double> coords;
                Vector1i1 tris;

                LoadObj3d(path.c_str(), coords, tris);

                if (coords.size() == 0)
                {
                    return;
                }

                for (int i = 0; i < (int)coords.size(); i += 3)
                {
                    vecs.push_back(Vector3d(coords[i + 0], coords[i + 1], coords[i + 2]));
                }

                for (int i = 0; i < (int)tris.size(); i += 3)
                {
                    face_id_0.push_back(tris[i + 0]);
                    face_id_1.push_back(tris[i + 1]);
                    face_id_2.push_back(tris[i + 2]);
                }
                /*********************************************************************************/
            }
        };

        //Output the vertex coordinates of a two-dimensional rectangle to a file in the specified path (one face)
        //One - dimensional vertex coordinates
        static void OutputRectangle2d(const std::string& path, const std::vector<Vector2d>& points)
        {
            std::ofstream file(path);

            for (int i = 0; i < points.size(); i++)
            {
                file << "v " << points[i][0] << " " << points[i][1] << " " << 0.0 << std::endl;//Output vertex coordinates(set z to zero)
            }

            int nb = 1;

            file << "f ";
            for (int i = 0; i < points.size(); i++)
            {
                file << Int2String(nb) << " ";//Output vertex number of face,number starts from 1
                nb++;
            }
            file << "" << std::endl;

            file.clear();
            file.close();
        }
        //Output the vertex coordinates  of 3D objects to a file in the specified path (one face)
        //One - dimensional vertex coordinates
        static void OutputObj3d(const std::string& path, const Vector3d1& points)
        {
            Functs::MAssert(points.size() >= 3,"CGAL_Output_Obj error: vecs.size() < 3 ");
            
            std::ofstream file(path);
            for (auto& p : points)
                file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;

            int nb = 1;
            file << "f ";
            for (int p = 0; p < points.size(); p++)
            {
                file << Int2String(nb) << " ";//Output vertex number of face
                nb++;
            }
            file << "" << std::endl;

            file.clear();
            file.close();
        };
        //Output the vertex coordinates and color values of 3D objects to a file in the specified path(one face)
        //One - dimensional vertex coordinates
        static void OutputObj3d(const std::string& path, const Vector3d1& points, const Vector3d1& colors)
        {
            Functs::MAssert(points.size()==colors.size(),"points.size()!=colors.size()");
            Functs::MAssert(points.size() >= 3,"CGAL_Output_Obj error: vecs.size() < 3 ");


            std::ofstream file(path);
            for(int i=0;i<points.size();i++)
            {
                auto p = points[i], c = colors[i];
                file << "v ";
                file << p[0] << " " << p[1] << " " << p[2] <<" ";
                file << c[0] << " " << c[1] << " " << c[2];
                file << std::endl;
            }
          
            int nb = 1;
            file << "f ";
            for (int p = 0; p < points.size(); p++)
            {
                file << Int2String(nb) << " ";
                nb++;
            }
            file << "" << std::endl;

            file.clear();
            file.close();
        };
        
        //Output the vertex coordinates  of 3D objects to a file in the specified path (multiple faces)
        //two - dimensional vertex coordinates
        static void OutputObj3d(const std::string& path, const Vector3d2& points, const int output_index = 1, const string str = "")
        {
            std::ofstream file(path);

            for (auto& points_ : points)
            {
                for (auto& p : points_)
                {
                    file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
                }
            }

            if (output_index == 2) file << "g " + str + "_p_0" << std::endl;

            int nb = 1;
            for (int i = 0; i < points.size(); i++)
            {
                if (output_index == 1)
                    file << "g " + str + "_f_" << i << std::endl;
                auto points_ = points[i];
                file << "f ";
                for (int p = 0; p < points_.size(); p++)
                {
                    file << Int2String(nb) << " ";
                    nb++;
                }
                file << "" << std::endl;
            }

            file.clear();
            file.close();
        };
        //Output the vertex coordinates  of 3D objects to a file in the specified path (multiple faces)
        //three - dimensional vertex coordinates
        static void OutputObj3d(const std::string& path, const Vector3d3& points, const int output_index = 1)
        {
            std::ofstream file(path);

            for (auto& points_ : points)
            {
                for (auto& points__ : points_)
                {
                    for (auto& p : points__)
                    {
                        file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
                    }
                }
            }

            int nb = 1;
            if (output_index == 3) file << "g ps_0" << std::endl;
            for (int i = 0; i < points.size(); i++)
            {
                if (output_index == 2) file << "g p_" << i << std::endl;

                for (int j = 0; j < points[i].size(); j++)
                {
                    if (output_index == 1) file << "g f_" << i << "_" << j << std::endl;

                    auto points_ = points[i][j];
                    file << "f ";
                    for (int p = 0; p < points_.size(); p++)
                    {
                        file << Int2String(nb) << " ";
                        nb++;
                    }
                    file << "" << std::endl;
                }
            }
            file.clear();
            file.close();
            
        };
        //Output the vertex coordinates and color values of 3D objects to a file in the specified path(multiple faces)
        //three - dimensional vertex coordinates,one-dimensional color value
        static void OutputObj3d(const std::string& path, const Vector3d3& points, const Vector3d1& colors, const int& output_index = 1)
        {
            std::ofstream file(path);

            for (int i = 0; i < points.size(); i++)
            {
                auto color = colors[i];
                for (int j = 0; j < points[i].size(); j++)
                {
                    for (int k = 0; k < points[i][j].size(); k++)
                    {
                        file << "v " << points[i][j][k][0] << " " << points[i][j][k][1] << " " << points[i][j][k][2] << " " << color[0] << " " << color[1] << " " << color[2] << std::endl;
                    }
                }
            }

            int nb = 1;
            if (output_index == 3) file << "g ps_0" << std::endl;
            for (int i = 0; i < points.size(); i++)
            {
                if (output_index == 2) file << "g p_" << i << std::endl;

                for (int j = 0; j < points[i].size(); j++)
                {
                    if (output_index == 1) file << "g f_" << i << "_" << j << std::endl;

                    auto points_ = points[i][j];
                    file << "f ";
                    for (int p = 0; p < points_.size(); p++)
                    {
                        file << Int2String(nb) << " ";
                        nb++;
                    }
                    file << "" << std::endl;
                }
            }
            file.clear();
            file.close();
        };
        //Output the vertex coordinates and color values of 3D objects to a file in the specified path(multiple faces)
        //three - dimensional vertex coordinates,two-dimensional color value
        static void OutputObj3d(const std::string& path, const Vector3d3& points, const Vector3d2& colors, int output_index = 1)
        {
            std::ofstream file(path);

            for (int i = 0; i < points.size(); i++)
            {
                for (int j = 0; j < points[i].size(); j++)
                {
                    auto color = colors[i][j];
                    for (int k = 0; k < points[i][j].size(); k++)
                    {
                        file << "v " << points[i][j][k][0] << " " << points[i][j][k][1] << " " << points[i][j][k][2] << " " << color[0] << " " << color[1] << " " << color[2] << std::endl;
                    }
                }
            }

            int nb = 1;
            if (output_index == 3) file << "g ps_0" << std::endl;
            for (int i = 0; i < points.size(); i++)
            {
                if (output_index == 2) file << "g p_" << i << std::endl;

                for (int j = 0; j < points[i].size(); j++)
                {
                    if (output_index == 1) file << "g f_" << i << "_" << j << std::endl;

                    auto points_ = points[i][j];
                    file << "f ";
                    for (int p = 0; p < points_.size(); p++)
                    {
                        file << Int2String(nb) << " ";
                        nb++;
                    }
                    file << "" << std::endl;
                }
            }
            file.clear();
            file.close();
        };
        //The parameter is LiblgpTriMesh
        static void OutputObj3d(const std::string& path, const LiblgpTriMesh& tm)
        {
            OutputObj3d(path, tm.vecs, tm.face_id_0, tm.face_id_1, tm.face_id_2);
        }
        //vecs:the vertex coordinates
        //face_id_0(1)(2): the index of triangular facesc:
        static void OutputObj3d(const std::string& path, const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2)
        {
            if (vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1)
            {
                std::cout << "vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1" << std::endl;
                return;
            }//Unable to form a triangular surface

            for (int i = 0; i < face_id_0.size(); i++)
            {
                int index_0 = face_id_0[i];
                int index_1 = face_id_1[i];
                int index_2 = face_id_2[i];
                if (index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size())
                {
                    std::cout << "index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size()" << std::endl;
                    return;
                }
            }//Vertex index out of bounds

            std::ofstream file(path);
            for (int i = 0; i < vecs.size(); i++)
            {
                OutputIterInfo("Output Vecs: ", (int)vecs.size(), i, 10);
                Vector3d v = vecs[i];
                file << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
            }

            for (int i = 0; i < face_id_0.size(); i++)
            {
                OutputIterInfo("Output faces: ", (int)face_id_0.size(), i, 10);
                int index_0 = face_id_0[i];
                int index_1 = face_id_1[i];
                int index_2 = face_id_2[i];
                if (index_0 != index_1 && index_0 != index_2 && index_1 != index_2)
                    file << "f " << index_0 + 1 << " " << index_1 + 1 << " " << index_2 + 1 << std::endl;
            }
            file.close();
        };

        static void OutputObj3d(const char* path, const Vector3d1& vecs, const std::vector<std::vector<int>>& face_ids)
        {
            std::vector<int> face_id_0;
            std::vector<int> face_id_1;
            std::vector<int> face_id_2;

            for (int i = 0; i < face_ids.size(); i++)
            {
                face_id_0.push_back(face_ids[i][0]);
                face_id_1.push_back(face_ids[i][1]);
                face_id_2.push_back(face_ids[i][2]);
            }

            OutputObj3d(path, vecs, face_id_0, face_id_1, face_id_2);
        }

        static void OutputObj3d(const char* path, const Vector3d1& vecs, const std::vector<std::vector<int>>& face_ids, const std::vector<int>& triangles_lables, const int& index)
        {
            std::vector<int> face_id_0;
            std::vector<int> face_id_1;
            std::vector<int> face_id_2;

            std::vector<int> lables(vecs.size(), -1);
            for (int i = 0; i < face_ids.size(); i++)
            {
                face_id_0.push_back(face_ids[i][0]);
                face_id_1.push_back(face_ids[i][1]);
                face_id_2.push_back(face_ids[i][2]);
                if (triangles_lables[i] == index)
                {
                    lables[face_ids[i][0]] = 0;
                    lables[face_ids[i][1]] = 0;
                    lables[face_ids[i][2]] = 0;
                }
            }
            Vector3d1 new_vecs;
            std::vector<int> new_face_id_0;
            std::vector<int> new_face_id_1;
            std::vector<int> new_face_id_2;

            int vertices_nb = 0;
            for (int i = 0; i < vecs.size(); i++)
            {
                if (lables[i] == 0)
                {
                    Vector3d v = vecs[i];
                    new_vecs.push_back(v);
                    lables[i] = vertices_nb;
                    vertices_nb++;
                }
            }

            for (int i = 0; i < face_id_0.size(); i++)
            {
                if (triangles_lables[i] == index)
                {
                    new_face_id_0.push_back(lables[face_id_0[i]]);
                    new_face_id_1.push_back(lables[face_id_1[i]]);
                    new_face_id_2.push_back(lables[face_id_2[i]]);
                }
            }

            OutputObj3d(path, new_vecs, new_face_id_0, new_face_id_1, new_face_id_2);
        }

        static void OutputObj3d(const char* path, const Vector3d1& vecs, const Vector3d1& colors, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2)
        {
            if (vecs.size() < 3 || colors.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1)
            {
                std::cout << "CGAL_Output_Obj error: vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1" << std::endl;
                return;
            }

            for (int i = 0; i < face_id_0.size(); i++)
            {
                int index_0 = face_id_0[i];
                int index_1 = face_id_1[i];
                int index_2 = face_id_2[i];

                if (index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size())
                {
                    std::cout << "CGAL_Output_Obj error: index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size()" << std::endl;
                    return;
                }
            }

            std::ofstream file(path);
            for (int i = 0; i < vecs.size(); i++)
            {
                file << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << " " << colors[i][0] << " " << colors[i][1] << " " << colors[i][2] << std::endl;
            }

            for (int i = 0; i < face_id_0.size(); i++)
            {
                int index_0 = face_id_0[i];
                int index_1 = face_id_1[i];
                int index_2 = face_id_2[i];

                if (index_0 != index_1 && index_0 != index_2 && index_1 != index_2)
                    file << "f " << index_0 + 1 << " " << index_1 + 1 << " " << index_2 + 1 << std::endl;
            }
            file.close();
        }

        static void OutputObj3d(const char* path, const Vector3d1& vecs, const Vector3d1& colors, const std::vector<std::vector<int>>& face_ids)
        {
            std::vector<int> face_id_0;
            std::vector<int> face_id_1;
            std::vector<int> face_id_2;
            for (int i = 0; i < face_ids.size(); i++)
            {
                face_id_0.push_back(face_ids[i][0]);
                face_id_1.push_back(face_ids[i][1]);
                face_id_2.push_back(face_ids[i][2]);
            }
            OutputObj3d(path, vecs, colors, face_id_0, face_id_1, face_id_2);
        }

        static void OutputOff3d(const char* path, const Vector3d1& vecs, const std::vector<int>& face_id_0, const std::vector<int>& face_id_1, const std::vector<int>& face_id_2)
        {
            if (vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1)
            {
                std::cout << "CGAL_Output_Off error: vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1" << std::endl;
                return;
            }

            for (int i = 0; i < face_id_0.size(); i++)
            {
                int index_0 = face_id_0[i];
                int index_1 = face_id_1[i];
                int index_2 = face_id_2[i];

                if (index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size())
                {
                    std::cout << "CGAL_Output_Off error: index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size()" << std::endl;
                    return;
                }
            }
            std::ofstream file(path);
            file << "OFF" << std::endl;
            file << vecs.size() << " " << face_id_0.size() << " 0" << std::endl;
            for (int i = 0; i < vecs.size(); i++)
                file << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
            for (int i = 0; i < face_id_0.size(); i++)
                file << "3 " << face_id_0[i] << " " << face_id_1[i] << " " << face_id_2[i] << " " << std::endl;
            file.close();
        }
        //Output color information to an. mtl file
        static void OutputMtl(const string& path, const Vector3d1& colors, const string& pre_name)
        {
            ofstream file(path);
            for (int i = 0; i < colors.size(); i++)
            {
                auto color = colors[i];
                //attribute line
                file << "newmtl " << pre_name << i << std::endl;
                file << "illum 4" << std::endl;
                file << "Kd " << color[0] << " " << color[1] << " " << color[2] << std::endl;
                file << "Ka 0.00 0.00 0.00" << std::endl;
                file << "Tf 1.00 1.00 1.00" << std::endl;
                file << "Ni 1.00" << std::endl;
            }
            file.close();
        }
        //Export geometric information to an output stream
        //e_index:Number of exported vertices
        //s_name:The name of the geometry
        //rgb:color value
        //mel_name:material quality

        static void ExportSegment(std::ofstream& output, int& e_index,
            const string& s_name, const Vector3d& start, const Vector3d& end,
            const double& radius, const bool cube_face = true)
        {
            ExportSegment(output, e_index, s_name, Vector3d(-1.0, -1.0, -1.0), "", start, end, radius, cube_face);
        }

        static void ExportSegment(std::ofstream& output, int& e_index,
            const string& s_name, const Vector3d& rgb, const Vector3d& start, const Vector3d& end,
            const double& radius, const bool cube_face = true)
        {
            ExportSegment(output, e_index, s_name, rgb, "", start, end, radius, cube_face);
        }

        static void ExportSegment(std::ofstream& output, int& e_index,
            const string& s_name, const string& mtl_name, const Vector3d& start, const Vector3d& end,
            const double& radius, const bool cube_face = true)
        {
            ExportSegment(output, e_index, s_name, Vector3d(0.0, 0.0, 0.0), mtl_name, start, end, radius, cube_face);
        }

        static void ExportSegment(std::ofstream& output, int& e_index,
            const string& s_name, const Vector3d& rgb, const string& mtl_name, const Vector3d& start, const Vector3d& end,
            const double& radius, const bool cube_face = true)
        {
            Vector3d normal = end - start;
            Vector3d base_1 = Vector3dBase(normal);

            base_1 = GetVectorLength(base_1, radius);

            Vector3d1 vecs;

            if (cube_face)
            {
                for (int i = 0; i < 4; i++) {
                    double angle = (double)(i)*Math_PI / 2.0;
                    Vector3d v = Functs::RotationAxis(normal + base_1, angle, normal);
                    vecs.push_back(v + start);
                }
                for (int i = 0; i < 4; i++) {
                    vecs.push_back(vecs[i] - normal);
                }
            }
            else
            {
                for (int i = 0; i < 2; i++) {
                    double angle = (double)(i)*Math_PI;
                    Vector3d v = Functs::RotationAxis(normal + base_1, angle, normal);
                    vecs.push_back(v + start);
                }
                for (int i = 0; i < 2; i++) {
                    vecs.push_back(vecs[i] - normal);
                }
            }

            Vector1i2 faces;
            if (cube_face)
            {
                int face_index_0[4] = { 0, 1, 2, 3 };
                int face_index_1[4] = { 5, 1, 0, 4 };
                int face_index_2[4] = { 4, 0, 3, 7 };
                int face_index_3[4] = { 5, 4, 7, 6 };
                int face_index_4[4] = { 7, 3, 2, 6 };
                int face_index_5[4] = { 6, 2, 1, 5 };

                faces.push_back(std::vector<int>(face_index_0, face_index_0 + 4));
                faces.push_back(std::vector<int>(face_index_1, face_index_1 + 4));
                faces.push_back(std::vector<int>(face_index_2, face_index_2 + 4));
                faces.push_back(std::vector<int>(face_index_3, face_index_3 + 4));
                faces.push_back(std::vector<int>(face_index_4, face_index_4 + 4));
                faces.push_back(std::vector<int>(face_index_5, face_index_5 + 4));

            }
            else
            {
                int face_index_0[4] = { 0, 1, 3, 2 };
                faces.push_back(std::vector<int>(face_index_0, face_index_0 + 4));
            }


            for (int i = 0; i < vecs.size(); i++)
            {
                if (mtl_name.size() != 0)
                {
                    output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;

                }
                else
                {
                    if (rgb[0] >= 0 && rgb[1] >= 0 && rgb[2] >= 0)
                        output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << " " <<
                        rgb[0] << " " << rgb[1] << " " << rgb[2] << std::endl;
                    else
                        output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
                }
            }

            if (std::string(s_name).size() != 0)
                output << "g " + std::string(s_name) << std::endl;

            if (mtl_name.size() != 0)
                output << "usemtl " << mtl_name << std::endl;

            for (int i = 0; i < faces.size(); i++) {
                output << "f ";
                for (int j = faces[i].size() - 1; j >= 0; j--)
                {
                    output << faces[i][j] + e_index << " ";
                }
                output << "" << std::endl;
            }

            e_index += (int)vecs.size();
        }
        //Export the geometric information of a stick to an output stream
        static void ExportStick(std::ofstream& output, int& e_index,
            const string& s_name, const Vector3d& rgb, const Vector3d1& start_poly, const Vector3d1& end_poly)
        {
            ExportStick(output, e_index, s_name, rgb, "", start_poly, end_poly);
        }

        static void ExportStick(std::ofstream& output, int& e_index,
            const string& s_name, const string& mtl_name, const Vector3d1& start_poly, const Vector3d1& end_poly)
        {
            ExportStick(output, e_index, s_name, Vector3d(-1.0, -1.0, -1.0), mtl_name, start_poly, end_poly);
        }

        static void ExportStick(std::ofstream& output, int& e_index,
            const string& s_name, const Vector3d1& start_poly, const Vector3d1& end_poly)
        {
            ExportStick(output, e_index, s_name, Vector3d(-1.0, -1.0, -1.0), "", start_poly, end_poly);
        }

        static void ExportStick(std::ofstream& output, int& e_index,
            const string& s_name, const Vector3d& rgb, const string& mtl_name, const Vector3d1& start_poly, const Vector3d1& end_poly)
        {
            //Functs::MAssert(start_poly.size()<3||end_poly.size()<3|| start_poly.size()!=end_poly.size(), "start_poly.size()<3||end_poly.size()<3|| start_poly.size()!=end_poly.size()");

            for (int i = 0; i < start_poly.size(); i++)
            {
                if (mtl_name.size() != 0)
                {
                    output << "v " << start_poly[i][0] << " " << start_poly[i][1] << " " << start_poly[i][2] << std::endl;

                }
                else
                {
                    if (rgb[0] >= 0 && rgb[1] >= 0 && rgb[2] >= 0)
                        output << "v " << start_poly[i][0] << " " << start_poly[i][1] << " " << start_poly[i][2] << " " <<
                        rgb[0] << " " << rgb[1] << " " << rgb[2] << std::endl;
                    else
                        output << "v " << start_poly[i][0] << " " << start_poly[i][1] << " " << start_poly[i][2] << std::endl;
                }
            }

            for (int i = 0; i < end_poly.size(); i++)
            {
                if (mtl_name.size() != 0)
                {
                    output << "v " << end_poly[i][0] << " " << end_poly[i][1] << " " << end_poly[i][2] << std::endl;

                }
                else
                {
                    if (rgb[0] >= 0 && rgb[1] >= 0 && rgb[2] >= 0)
                        output << "v " << end_poly[i][0] << " " << end_poly[i][1] << " " << end_poly[i][2] << " " <<
                        rgb[0] << " " << rgb[1] << " " << rgb[2] << std::endl;
                    else
                        output << "v " << end_poly[i][0] << " " << end_poly[i][1] << " " << end_poly[i][2] << std::endl;
                }
            }


            if (std::string(s_name).size() != 0)
                output << "g " + std::string(s_name) << std::endl;

            if (mtl_name.size() != 0)
                output << "usemtl " << mtl_name << std::endl;

            output << "f ";
            for (int i = 0; i < start_poly.size(); i++)
                output << i + e_index << " ";
            output << std::endl;

            output << "f ";
            for (int i = end_poly.size() - 1; i >= 0; i--)
                output << i + e_index + start_poly.size() << " ";
            output << std::endl;

            for (int i = 0; i < start_poly.size(); i++)
            {
                auto a = std::to_string(e_index + i);
                auto b = std::to_string(e_index + (i + 1) % start_poly.size());
                auto a_ = std::to_string(start_poly.size() + e_index + i);
                auto b_ = std::to_string(start_poly.size() + e_index + (i + 1) % start_poly.size());

                output << "f " << a_ << " " << b_ << " " << b << " " << a << std::endl;;
            }

            e_index += (int)start_poly.size() + (int)end_poly.size();
        }
        //Export the vertex coordinates and face indices of a cube out of the output file stream to represent a point
        static void ExportPoint(std::ofstream& output, int& e_index, const string& s_name, const Vector3d point, const double& radius)
        {
            ExportPoint(output, e_index, s_name, Vector3d(-1.0, -1.0, -1.0), "", point, radius);
        }

        static void ExportPoint(std::ofstream& output, int& e_index, const string& s_name, const Vector3d& rgb, const Vector3d point, const double& radius)
        {
            ExportPoint(output, e_index, s_name, rgb, "", point, radius);
        }

        static void ExportPoint(std::ofstream& output, int& e_index, const string& s_name, const string& mtl_name, const Vector3d point, const double& radius)
        {
            ExportPoint(output, e_index, s_name, Vector3d(0.0, 0.0, 0.0), mtl_name, point, radius);
        }

        static void ExportPoint(std::ofstream& output, int& e_index, const Vector3d point, const double& radius)
        {
            ExportPoint(output, e_index, "", Vector3d(-1.0, -1.0, -1.0), "", point, radius);
        }

        static void ExportPoint(std::ofstream& output, int& e_index, const Vector3d& rgb, const Vector3d point, const double& radius)
        {
            ExportPoint(output, e_index, "", rgb, "", point, radius);
        }
        
        static void ExportPoint
        (std::ofstream& output, int& e_index, const string& s_name,
            const Vector3d& rgb, const string& mtl_name, const Vector3d point, const double& radius)
        {
            Vector3d1 vecs;//The set of vertices of a unit cube
            vecs.push_back(Vector3d(0.5, 0.5, 0.5));
            vecs.push_back(Vector3d(-0.5, 0.5, 0.5));
            vecs.push_back(Vector3d(-0.5, 0.5, -0.5));
            vecs.push_back(Vector3d(0.5, 0.5, -0.5));

            vecs.push_back(Vector3d(0.5, -0.5, 0.5));
            vecs.push_back(Vector3d(-0.5, -0.5, 0.5));
            vecs.push_back(Vector3d(-0.5, -0.5, -0.5));
            vecs.push_back(Vector3d(0.5, -0.5, -0.5));

            Vector1i2 faces;

            int face_index_0[4] = { 0, 1, 2, 3 };
            int face_index_1[4] = { 5, 1, 0, 4 };
            int face_index_2[4] = { 4, 0, 3, 7 };
            int face_index_3[4] = { 5, 4, 7, 6 };
            int face_index_4[4] = { 7, 3, 2, 6 };
            int face_index_5[4] = { 6, 2, 1, 5 };

            faces.push_back(std::vector<int>(face_index_0, face_index_0 + 4));
            faces.push_back(std::vector<int>(face_index_1, face_index_1 + 4));
            faces.push_back(std::vector<int>(face_index_2, face_index_2 + 4));
            faces.push_back(std::vector<int>(face_index_3, face_index_3 + 4));
            faces.push_back(std::vector<int>(face_index_4, face_index_4 + 4));
            faces.push_back(std::vector<int>(face_index_5, face_index_5 + 4));

            for (int i = 0; i < vecs.size(); i++) {
                vecs[i][0] = vecs[i][0] * radius;//Scale
                vecs[i][1] = vecs[i][1] * radius;
                vecs[i][2] = vecs[i][2] * radius;

                vecs[i][0] += point[0];//translation
                vecs[i][1] += point[1];
                vecs[i][2] += point[2];

                if (mtl_name.size() != 0)
                    output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
                else
                {
                    if (rgb[0] >= 0 && rgb[1] >= 0 && rgb[2] >= 0)
                        output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << " " << rgb[0] << " " << rgb[1] << " " << rgb[2] << std::endl;
                    else
                        output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
                }
            }

            if (std::string(s_name).size() != 0)
                output << "g " + std::string(s_name) << std::endl;

            if (mtl_name.size() != 0)
                output << "usemtl " << mtl_name << std::endl;


            for (int i = 0; i < faces.size(); i++) {
                output << "f ";

                for (int j = (int)faces[i].size() - 1; j >= 0; j--) {
                    output << faces[i][j] + e_index << " ";
                }
                output << "" << std::endl;
            }
            e_index += 8;
        }
        //Output the information of the tree to a file at the given path
        static void Outputtree(const int& nodes_nb, const std::vector<int>& edges,
            const std::string& path, const std::vector<string> labels = std::vector<string>())
        {
            std::ofstream file(path);

            file << "Mark Newman on Sat Jul 22 05:32:16 2006" << std::endl;
            file << "graph" << std::endl;
            file << "[" << std::endl;
            file << "  directed 0" << std::endl;

            for (int i = 0; i < nodes_nb; i++)
            {
                file << "node" << std::endl;
                file << "[" << std::endl;
                file << "id " << i << std::endl;

                if (labels.size() == nodes_nb)
                    file << "label " << labels[i] << std::endl;
                else
                    file << "label " << i << std::endl;

                file << "]" << std::endl;
            }

            for (int i = 0; i < edges.size(); i = i + 2)
            {
                file << "edge" << std::endl;
                file << "[" << std::endl;

                file << "source " << edges[i] << std::endl;
                file << "target " << edges[static_cast<int64_t>(i) + 1] << std::endl;

                file << "]" << std::endl;
            }

            file << "]" << std::endl;

            file.clear();
            file.close();
        }


#ifndef __APPLE__
        //Clear all files and folders in the specified path
        //If the path does not exist, the function will create all folders in the path
        //If the path exists, first delete all files and folders under the path, and then create a new empty folder
        static void ClearFolder(const std::string& path)
        {
			if (!DetectExisting(path))
			{
				auto folders = SplitStr(StringReplace(path, "\\", "/"), "/");
				if (folders.back() == "")
					folders.erase(folders.begin() + folders.size() - 1);
				std::string str;
				for (int i = 0; i < folders.size(); i++)
				{
					str += folders[i] + "/";
					if (!DetectExisting(str))
					{
						if (_mkdir(str.c_str())) {};//create the floder
					}
				}
			}
			else
			{
				std::string del_cmd = "del /f/s/q " + path + " > nul";
				system(del_cmd.c_str());//Delete all files under the path
				std::string rmdir_cmd = "rmdir /s/q " + path;
				system(rmdir_cmd.c_str());//Delete all folders under the path

				if (_mkdir(path.c_str())) {};
			}
        }
        //Detect if the specified path exists
        static bool DetectExisting(const std::string& path)
        {
            struct stat buffer;
            return (stat(path.c_str(), &buffer) == 0);
        }
#else
        static void ClearFolder(const std::string& path)
        {
            if (access(path.c_str(), 0) == -1)
            {
				auto folders = SplitStr(StringReplace(path, "\\", "/"), "/");
				if (folders.back() == "")
					folders.erase(folders.begin() + folders.size() - 1);
				std::string str;
				for (int i = 0; i < folders.size(); i++)
				{
					str += folders[i] + "/";
					if (!DetectExisting(str))
					{
                        system(std::string("mkdir " + str).c_str());
					}
				}
            }
            else
            {
				std::string del_cmd = "rm -rf " + path;
				system(del_cmd.c_str());
				std::string mkdir_cmd = "mkdir " + path;
				system(mkdir_cmd.c_str());
            }
        }
        
        static bool DetectExisting(const std::string& path)
        {
            if (access(path.c_str(), 0) == -1)
                return false;
            else
                return true;
        }
#endif



        //Generate a complete file path based on the executable file path of the program and the parameters passed in and detect its existence
#ifndef __APPLE__
        static std::string GenerateFullPath(const std::string& py_path)
        {
            std::string path = std::string(_pgmptr).substr(0, std::string(_pgmptr).find_last_of('\\')) + py_path;
            //std::string path = std::string(_pgmptr).substr(0, std::string(_pgmptr).find_last_of('\\')) + py_path;
            if (!Functs::DetectExisting(path))
                Functs::MAssert("std::string EXP(const std::string py_path="")");
            return path;
        }
        //Obtain and return all file names in the specified directory
        static VectorStr1 GetFilesInDirectory(const std::string& path)
        {

            vector<string> names;

#if !defined(UNICODE) && !defined(_UNICODE)
            string search_path = path + "/*.*";
            WIN32_FIND_DATA fd;
            HANDLE hFind = ::FindFirstFile(search_path.c_str(), &fd);
            if (hFind != INVALID_HANDLE_VALUE) {
                do {
                    // read all (real) files in current folder
                    // , delete '!' read other 2 default folder . and ..
                    if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
                        names.push_back(fd.cFileName);
                    }
                } while (::FindNextFile(hFind, &fd));
                ::FindClose(hFind);
            }
#endif

            return names;
        }
#else
        static VectorStr1 GetFilesInDirectory(const std::string& path)
        {
            vector<string> names;
            //for (const auto & entry : fs::directory_iterator(path))
             //   names.push_back(entry.path());
            Functs::MAssert("static VectorStr1 GetFilesInDirectory(const std::string& path)");
            return names;
        }
#endif



#pragma endregion

#pragma region Graph
        //Calculate the minimum spanning tree
        //edges:Elements are nodes, and adjacent elements represent the starting and ending points of an edge
        static std::vector<int> MinimalSpanningTreeGeneral(
            const std::vector<int>& edges,
            const std::vector<double>& costs)
        {
            auto unique_nodes = UniqueSet(edges);

            int node_nb = static_cast<int>(unique_nodes.size());

            std::map<int, int> unique_map_0;
            std::map<int, int> unique_map_1;
            for (int i = 0; i < unique_nodes.size(); i++)
            {
                unique_map_0.insert(std::pair<int, int>(unique_nodes[i], i));
                unique_map_1.insert(std::pair<int, int>(i, unique_nodes[i]));
            }

            std::vector<int> map_edges;
            for (auto edge : edges) map_edges.emplace_back(unique_map_0.at(edge));//Map nodes in edges to indexes

            auto mst = MinimalSpanningTree(node_nb, map_edges, costs);

            std::vector<int> map_mst;
            for (auto m : mst) map_mst.emplace_back(unique_map_1.at(m));//Remap index to node

            return map_mst;
        }
        //Calculate the minimum spanning tree
        //edges:The element is the index of a node, and adjacent elements represent the starting and ending points of an edge
        static std::vector<int> MinimalSpanningTree(
            const int& node_nb, const std::vector<int>& edges, const std::vector<double>& costs)
        {
            std::vector<int> mst;
            std::vector<int> nodes;

            for (int i = 0; i < node_nb; i++) nodes.push_back(i);

            std::vector<std::vector<int>> containers;
            for (int i = 0; i < nodes.size(); i++)
            {
                std::vector<int> container;
                container.push_back(nodes[i]);
                containers.push_back(container);
                std::vector<int>().swap(container);//Free up empty space
            }

            std::vector<bool> edges_used;
            for (int i = 0; i < costs.size(); i++)
                edges_used.push_back(false);//Usage of storage edges

            do
            {
                //find a minimal cost edge
                int minimal_cost_edge_index = -1;
                double minimal_cost = MAXDOUBLE;
#pragma region find_a_minimal_cost_edge

                for (int j = 0; j < costs.size(); j++)
                {
                    if (!edges_used[j])
                    {
                        if (costs[j] < minimal_cost)
                        {
                            minimal_cost = costs[j];
                            minimal_cost_edge_index = j;
                        }
                    }
                }
#pragma endregion

                if (minimal_cost_edge_index < 0)
                    break;

                //check valid
                int edge_index_0 = static_cast<int>(2)* minimal_cost_edge_index;//Calculate the index of the starting node of the minimum cost edge in edges
                int edge_index_1 = static_cast<int>(2)* minimal_cost_edge_index + 1; //Calculate the index of the termination node of the minimum cost edge in edges
                int node_index_0 = edges[edge_index_0];
                int node_index_1 = edges[edge_index_1];

                int container_0 = -1;
                int container_0_0 = -1;
                int container_1 = -1;
                int container_1_0 = -1;

                for (int j = 0; j < containers.size() && (container_0 < 0 || container_1 < 0); j++)
                {
                    for (int k = 0; k < containers[j].size() && (container_0 < 0 || container_1 < 0); k++)
                    {
                        if (node_index_0 == containers[j][k])
                        {
                            container_0 = j;
                            container_0_0 = k;
                        }
                        if (node_index_1 == containers[j][k])
                        {
                            container_1 = j;
                            container_1_0 = k;
                        }
                    }
                }

                if (!(container_0 >= 0 && container_1 >= 0))
                {
                    break;
                }

                if (container_0 == container_1)
                {
                    edges_used[minimal_cost_edge_index] = true;
                }
                else
                {
                    mst.push_back(node_index_0);
                    mst.push_back(node_index_1);
                    edges_used[minimal_cost_edge_index] = true;

                    for (int i = 0; i < containers[container_1].size(); i++)
                    {
                        containers[container_0].push_back(containers[container_1][i]);
                    }

                    containers.erase(containers.begin() + container_1);
                }

            } while (containers.size() != 1);

            std::vector<bool>().swap(edges_used);
            std::vector<std::vector<int>>().swap(containers);

            std::vector<int>().swap(nodes);


            return mst;
        };
        //Find all connected components of a given tree
        //tree:The element is the index of the starting and ending nodes of each edge
        //Returns the index of the node in the connected component
        static std::vector<std::vector<int>> ConnectedComponents(const int& node_nb, const std::vector<std::pair<int, int>>& tree)
        {
            std::vector<int> tree_;
            for (auto& o : tree)
            {
                tree_.emplace_back(o.first);
                tree_.emplace_back(o.second);
            }
            return ConnectedComponents(node_nb, tree_);
        }
        //Find all connected components of a given tree
        //tree:Store node index with an edge between adjacent nodes
        //Returns the index of the node in the connected component
        static std::vector<std::vector<int>> ConnectedComponents(const int& node_nb, const std::vector<int>& tree)
        {
            std::vector<std::vector<int>> components;//Store all connected components
            std::vector<int> index(node_nb, -1);//Index storing the connected components to which each node belongs
            int nb = 0;
            for (int i = 0; i < tree.size(); i = i + 2)
            {
                int ii = i + 1;
                if (index[tree[i]] == -1 && index[tree[ii]] == -1)
                {
                    index[tree[i]] = nb;
                    index[tree[ii]] = nb;
                    nb++;
                }
                if (index[tree[i]] == -1 && index[tree[ii]] != -1)
                {
                    index[tree[i]] = index[tree[ii]];
                }
                if (index[tree[i]] != -1 && index[tree[ii]] == -1)
                {
                    index[tree[ii]] = index[tree[i]];
                }

                if (index[tree[i]] != -1 && index[tree[ii]] != -1)
                {
                    int min_index = std::min(index[tree[i]], index[tree[ii]]);
                    int max_index = std::max(index[tree[i]], index[tree[ii]]);
                    for (auto& index_ : index)
                    {
                        if (index_ == max_index)index_ = min_index;
                    }
                }

            }

            for (int i = 0; i < nb; i++)
            {
                std::vector<int> one;
                for (int j = 0; j < index.size(); j++)
                    if (index[j] == i)
                        one.emplace_back(j);
                if (!one.empty())components.emplace_back(one);
            }

            for (int j = 0; j < index.size(); j++) if (index[j] == -1)components.emplace_back(std::vector<int>(1, j));//Add nodes that are not assigned to connected components as a separate connected component to components

            for (auto& component : components)
                std::sort(component.begin(), component.end());
            return components;
        };
        //Find all connected components of a given tree
        //Returns the node in the connected component
        //tree:The element is the index of the starting and ending nodes of each edge
        static std::vector<std::vector<int>> ConnectedComponentsGeneral(const Vector1i1& nodes, const std::vector<std::pair<int, int>>& tree)
        {
            std::vector<int> tree_;
            for (auto& o : tree)
            {
                tree_.emplace_back(o.first);
                tree_.emplace_back(o.second);
            }
            return ConnectedComponentsGeneral(nodes, tree_);
        }
        //tree:Store node index with an edge between adjacent nodes
        static std::vector<std::vector<int>> ConnectedComponentsGeneral(const Vector1i1& nodes, const std::vector<int>& tree)
        {
            int node_nb = static_cast<int>(nodes.size());

            std::map<int, int> unique_map_0;
            std::map<int, int> unique_map_1;
            for (int i = 0; i < nodes.size(); i++)
            {
                unique_map_0.insert(std::pair<int, int>(nodes[i], i));
                unique_map_1.insert(std::pair<int, int>(i, nodes[i]));
            }

            std::vector<int> map_edges;
            for (auto edge : tree) map_edges.emplace_back(unique_map_0.at(edge));//Map nodes in edges to indexes

            auto components = ConnectedComponents(node_nb, map_edges);

            std::vector<std::vector<int>> map_components;
            for (auto component : components)
            {
                map_components.emplace_back(std::vector<int>());
                for (auto c : component)
                    map_components.back().emplace_back(unique_map_1.at(c));//Remap index to node
            }
            return map_components;
        }

        

#pragma endregion

#pragma region DevelopmentRelated
        //Print a line of text in the standard error stream (std:: cerr)
        //Returns true if printing is successful
        static bool CerrLine(const string& line, const int level = 0)
        {
            for (int i = 0; i < level; i++)
                std::cerr << CERR_ITER;
            std::cerr << line << std::endl;
            return true;
        }
        //Print a line of text in an output file stream (ofstream) and standard error stream (std:: cerr)
        //Returns true if printing is successful
        static bool CerrLine(ofstream& file, const std::string& line, const int level = 0)
        {
            for (int i = 0; i < level; i++)
            {
                file << CERR_ITER;
                std::cerr << CERR_ITER;
            }
            file << line << std::endl;
            std::cerr << line << std::endl;
            return true;
        }

        //Output iteration information in standard error stream (std:: cerr)
        //tn: total number of iterations
        //cn: current iteration
        //fn: output frequency number
        static void OutputIterInfo(const string& title, const int& tn, const int& cn, const int& fn, const int level = 0)
        {
            int delta = fn > tn ? 1 : tn / fn;

            if (cn % delta == 0)
            {
                if (cn == 0)
                {
                    for (int i = 0; i < level; i++)
                        std::cerr << CERR_ITER;
                    std::cerr << title << ": ";
                }
                std::cerr << Functs::Double2String((double)(100.0 * cn / tn), 1) << "% ";
                if (cn + delta >= tn) std::cerr << std::endl;
            }
        }
        //Output an error message in the standard error stream (std:: cerr) and perform some actions as needed
        //MacOS:using the 'read - p' command,press Enter to continuepress any key to continue
        //Windows:using the 'pause' command,press any key to continue
        static void MAssert(const std::string& str, const double sleep_seconds = -1)
        {
            std::cerr << "Bug: " << str << std::endl;
            if (sleep_seconds <= 0)
            {
                if (MACEN)system("read -p 'Press Enter to continue...' var");
                if (WINEN)system("pause");
            }
            else
                MSleep(sleep_seconds);
        }

        static void MAssert(const char* str, const double sleep_seconds = -1)
        {
            MAssert(std::string(str), sleep_seconds);
        }


        static bool MAssert(const bool& b, const std::string& str, const double sleep_seconds = -1)
        {
            if (!b)MAssert(str, sleep_seconds);
            return b;
        }

        static bool MAssert(const bool& b, const char* str, const double sleep_seconds = -1)
        {
            if (!b)MAssert(str, sleep_seconds);
            return b;
        }
        //Let the current thread sleep for a period of time, in seconds
        static void MSleep(const double& second)
        {
            this_thread::sleep_for(chrono::milliseconds((int)(second * 1000)));
        }


#ifndef __APPLE__
        //Running Python scripts in C++
        static void RunPY(const std::string& py_path, const std::string& paras)
        {
            std::string cmd = "python " + Functs::GenerateFullPath(py_path) + " " + paras;
            system(cmd.c_str());
        }
        //Obtain the path of the current working directory and return it as a string
        static std::string WinGetCurDirectory()
        {
            char tmp[256];
            if (_getcwd(tmp, 256)) {};
            return std::string(tmp);
        }
#else
        static std::string WinGetCurDirectory()
        {
            char tmp[256];
            if (getcwd(tmp, 256)) {};
            return std::string(tmp);
        }
#endif
        //Running system commands in C++
        static void RunCMD(const std::string& cmd_str)
        {
            std::cerr << "Command String: " << cmd_str << std::endl;;
            system(cmd_str.c_str());
        }
        //Obtain the username of the current user and return it as a string
        static std::string WinGetUserName()
        {
            char* user = getenv("username");
            return std::string(user);
        }

        //Copying files on Windows or Mac systems
        //source_file:Source file path
        //target_floder:target_folder path
        static bool WinCopy(const std::string& source_file, const std::string& target_folder, const bool& b = true)
        {
            if (!Functs::DetectExisting(source_file))
            {
                MAssert("Source file does not exist: " + source_file);
                return false;
            }

            if (!Functs::DetectExisting(target_folder))
            {
                MAssert("Target folder does not exist: " + target_folder);
                return false;
            }

            if (WINEN)
            {
                std::string str = "copy " + source_file + " " + target_folder;
                if (b) std::cerr << "Command string: " << str << std::endl;
                system(str.c_str());
            }
            if (MACEN)
            {
                //% cp -R ~/Documents/Expenses /Volumes/Data/Expenses
                std::string str = "cp -R " + source_file + " " + target_folder;
                if (b) std::cerr << "Command string: " << str << std::endl;
                system(str.c_str());
            }
            return true;
        }
        //Delete files on Windows system
        static bool WinDel(const std::string& source_file, const bool& b = true)
        {
            if (!Functs::DetectExisting(source_file))
            {
                //MAssert("Source file does not exist: " + source_file);
                return false;
            }
            std::string str = "del " + source_file;
            if (b) std::cerr << "Command string: " << str << std::endl;
            system(str.c_str());
            return true;
        }
        //Rename files in Windows systems
        static bool WinRename(const std::string& source_file, const std::string& rename_file, const bool& b = true)
        {
            if (!Functs::DetectExisting(source_file))
            {
                MAssert("Source file does not exist: " + source_file);
                return false;
            }

            std::string str = "rename " + source_file + " " + rename_file;
            if (b) std::cerr << "Command string: " << str << std::endl;
            system(str.c_str());
            return true;
        }
        //Obtain the type information of the parameter and return it as a string
        template <class Type>
        static std::string GetTypeId(const Type& t)
        {
            return  typeid(t).name();
        }

#pragma endregion
        //Returns a three-dimensional vector representing color(RGB)
        static Vector3d ColorMapping(const int& cur, const int& all)
        {
            MAssert(cur >= 0 && all >= 0 && cur <= all - 1, "cur>=0 && all>=0 && cur<all");
            double isolevel = (double)cur / (double)all;
            return ColorMapping(isolevel);
        }

        static Vector3d ColorMapping(const double& isolevel)
        {
            double output_c_0, output_c_1, output_c_2;
            ColorMapping(isolevel, output_c_0, output_c_1, output_c_2);
            return Vector3d(output_c_0, output_c_1, output_c_2);
        }

        static void ColorMapping(const double& isolevel, double& output_c_0, double& output_c_1, double& output_c_2)
        {
            MAssert(isolevel >= 0.0 && isolevel <= 1.0, "isolevel>=0.0&&isolevel<=1.0");

            Vector3d v;
            if (isolevel >= 0 && isolevel <= 0.25)
            {
                v[0] = 0;
                v[1] = isolevel / 0.25;
                v[2] = 1;
            }

            if (isolevel > 0.25 && isolevel <= 0.50)
            {
                v[0] = 0;
                v[1] = 1;
                v[2] = 1 - (isolevel - 0.25) / 0.25;
            }

            if (isolevel > 0.50 && isolevel <= 0.75)
            {
                v[0] = (isolevel - 0.50) / 0.25;
                v[1] = 1;
                v[2] = 0;
            }

            if (isolevel > 0.75 && isolevel <= 1.0)
            {
                v[0] = 1;
                v[1] = 1 - (isolevel - 0.75) / 0.25;
                v[2] = 0;
            }

            if (isolevel < 0.0)
            {
                v[0] = 0.0;
                v[1] = 0.0;
                v[2] = 0.0;
            }

            if (isolevel > 1.0)
            {
                v[0] = 0.5;
                v[1] = 0.0;
                v[2] = 0.0;
            }
            output_c_0 = v[0];
            output_c_1 = v[1];
            output_c_2 = v[2];
        }
    };

    typedef Functs FF;
    #define DS FF::DoubleString
    #define IS FF::IntString
	
    //To debug a release build
    //Open the Property Pages dialog box for the project.
    //Click the C / C++ node. Set Debug Information Format to C7 compatible(/ Z7) or Program Database(/ Zi).
    //Expand Linker and click the General node.Set Enable Incremental Linking to No(/ INCREMENTAL:NO).
    //Select the Debugging node.Set Generate Debug Info to Yes(/ DEBUG).
    //Select the Optimization node.Set References to / OPT:REF and Enable COMDAT Folding to / OPT : ICF.


}
#endif
