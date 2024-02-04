#ifndef _RI_H
#define _RI_H

#include <iostream>
#include <map>
#include "liblgp.hpp"
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>

using namespace liblgp;


//Resume IO Info
struct RI
{
public:

    int numb = 0;
    std::ofstream ofile;
    std::ifstream ifile;
    bool ol;//TRUE:The file path does not exist£»FALSE:The file path  exists
    bool cout_log=false;
    //Constructor
    RI(const string path) :numb(0), cout_log(false)
    {
        ol = !Functs::DetectExisting(path);

        if (ol)//The file path does not exist
        {
            ofile.open(path);//Automatically create files and open them
        }
        else
        {
            if (!Functs::DetectExisting(path))
                Functs::MAssert("void Init(const string path)");
            ifile.open(path);//If the file path exists, open the file directly
        }
    }


    RI(const bool ol_, const string path) :numb(0), ol(ol_), cout_log(false)
    {
        if (ol_)
        {
            ofile.open(path);
        }
        else
        {
            if (!Functs::DetectExisting(path))
                Functs::MAssert("void Init(const string path)");
            ifile.open(path);
        }
    }


    RI(const bool ol_, const string path, const bool cout_log_) :numb(0), ol(ol_), cout_log(cout_log_)
    {
        if (ol_)
        {
            ofile.open(path);
        }
        else
        {
            if (!Functs::DetectExisting(path))
                Functs::MAssert("void Init(const string path)");
            ifile.open(path);
        }
    }
    //close file
    void Close()
    {
        if (ol)
        {
            ofile.close();
            ofile.clear();
        }
        else
        {
            ifile.close();
            ifile.clear();
        }

        numb = 0;
    }

    int Check()
    {
        int nb = numb;
        if (nb == 703)
        {
            //Functs::MAssert("int Check()");
        }
        return nb;
    }

    void Read(const string t, const string l, const string& name)
    {
        if (l != name)
        {
            Functs::MAssert("Error: if (l != name)");
        }

        if (true)
        {
            if (t == "705:")
            {
                //int dsad = 0;
                //Functs::MAssert("577");
            }

            if (cout_log && atoi(t.c_str()) % 100 == 0)
                std::cerr << t << std::endl;
        }
    }

    int NB(const string& name, const int t)
    {
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t << std::endl;
            return t;
        }
        else
        {
            string tmp, lmp; int nb;
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            return nb;
        }
    }

    //Out:output numb,name and t,or read information from the file(depend on ol)

    //When ol is true (file does not exist), output "numb:name t" to the created file
    //When ol is false (file exists), read information from the file (stored in imp, tmp)
    void Out(const string& name, int& t)
    {
        string tmp, lmp;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t << std::endl;
        }
        else
        {
            ifile >> lmp >> tmp >> t;
            Read(lmp, tmp, name);
        }
    }
    
    void Out(const string& name, double& t)
    {
        string tmp, lmp;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();

            if (isnan(t))
            {
                ofile << name << " " << 1.7976931348623157E+308 << std::endl;
            }
            else
            {
                //std::numeric_limits<double>::infinity();
                if (t == std::numeric_limits<double>::infinity())
                {
                    //ofile << name << " " << std::numeric_limits<double>::max() << std::endl;
                    ofile << name << " " << 1.7976931348623157E+308 << std::endl;
                }
                else
                    ofile << fixed << setprecision(8) << name << " " << t << std::endl;
            }
        }
        else
        {
            ifile >> lmp >> tmp >> t;
            Read(lmp, tmp, name);
        }
    }

    void Out(const string& name, bool& t)
    {
        string tmp, lmp;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t << std::endl;
        }
        else
        {
            ifile >> lmp >> tmp >> t;
            Read(lmp, tmp, name);
        }
    }
   
    void Out(const string& name, string& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            if (t.size() == 0) ofile << name << " " << t.size() << std::endl;
            if (t.size() != 0) ofile << name << " " << t.size() << " " << t << std::endl;
        }
        else
        {
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            if (nb != 0) ifile >> t;
        }
    }

    //output path
    //pt:pref path
    void Out(const string& name, const string& pt, string& t)
    {
        if (ol)
        {
            if (t.size() == 0)
            {
                Out(name, t);
            }
            else
            {
                if (t.find(pt) == std::string::npos)//t does not include pt
                {
                    Functs::MAssert("if (t.find(pt) == std::string::npos)");
                }
                auto local_path = t.substr(pt.size(), t.size());
                Out(name, local_path);
            }
        }
        else
        {
            Out(name, t);
            if (t.size() != 0)
                t = pt + t;
        }
    }
  
    void Out(const string& name, VectorTI3& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << std::get<0>(t_) << " " << std::get<1>(t_) << " " << std::get<2>(t_) << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = VectorTI3(nb, std::tuple<int, int, int>());
            for (auto& t_ : t)
                ifile >> std::get<0>(t_) >> std::get<1>(t_) >> std::get<2>(t_);
        }
    }
   
    void Out(const string& name, Vector1i1& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++;
            ofile << numb << ": ";
            Check();
            ofile << name << " " << t.size() << " " << Functs::Int2String(t, false, " ") << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = Vector1i1(nb, -1);
            for (auto& t_ : t) ifile >> t_;
        }
    }
    
    void Out(const string& name, Vector1d1& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " " << Functs::Double2String(t, 8, false, " ") << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = Vector1d1(nb, -1);
            for (auto& t_ : t) ifile >> t_;
        }
    }
   
    void Out(const string& name, Vector1i2& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.size() << " " << Functs::Int2String(t_, false, " ") << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = Vector1i2(nb, Vector1i1());
            for (auto& t_ : t)
            {
                ifile >> nb;
                t_ = Vector1i1(nb, -1);
                for (auto& t__ : t_)
                    ifile >> t__;
            }
        }
    }

    void Out(const string& name, Vector1i3& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.size() << " ";
                for (auto& t__ : t_)
                {
                    ofile << t__.size() << " " << Functs::Int2String(t__, false, " ") << " ";
                }
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = Vector1i3(nb, Vector1i2());
            for (auto& t_ : t)
            {
                ifile >> nb;
                t_ = Vector1i2(nb, Vector1i1());
                for (auto& t__ : t_)
                {
                    ifile >> nb;
                    t__ = Vector1i1(nb,-1);
                    for(auto&t___:t__)
                        ifile >> t___;

                }
            }
        }
    }

    void Out(const string& name, Vector1d2& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.size() << " " << Functs::Double2String(t_, 8, false, " ") << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = Vector1d2(nb, Vector1d1());
            for (auto& t_ : t)
            {
                ifile >> nb;
                t_ = Vector1d1(nb, -1);
                for (auto& t__ : t_)
                    ifile >> t__;
            }
        }
    }
    void Out(const string& name, VectorPI1& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " " << Functs::Int2String(t, false, " ", " ") << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = VectorPI1(nb, std::pair<int, int>());
            for (auto& t_ : t)
                ifile >> t_.first >> t_.second;
        }
    }

    void Out(const string& name, const string& pt, std::vector<string>& t)
    {
        if (ol)
        {
            std::vector<string> local_paths;
            for (auto& t_ : t)
            {
                auto local_path = t_.substr(pt.size(), t_.size());
                local_paths.push_back(local_path);
            }
            Out(name, local_paths);
        }
        else
        {
            t.clear();
            std::vector<string> local_paths;
            Out(name, local_paths);
            for (auto& lp : local_paths)
            {
                t.push_back(pt + lp);
            }
        }
    }

    void Out(const string& name, std::vector<string>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.size() << " ";
                if (t_.size() != 0) ofile << t_ << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = std::vector<string>(nb, string(""));
            for (auto& t_ : t)
            {
                ifile >> nb;
                if (nb != 0) ifile >> t_;
            }
        }
    }
    void Out(const string& name, std::vector <std::vector<string>>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << std::endl;
            for (auto& t_ : t)
                Out(name, t_);
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                std::vector<string> t_;
                Out(name, t_);
                t.push_back(t_);
            }
        }
    }
    void Out(const string& name, VectorPI2& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
                ofile << t_.size() << " " << Functs::Int2String(t_, false, " ", " ") << " ";
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = VectorPI2(nb, VectorPI1());
            for (auto& t_ : t)
            {
                ifile >> nb;
                t_ = VectorPI1(nb, std::pair<int, int>());
                for (auto& t__ : t_)
                    ifile >> t__.first >> t__.second;
            }
        }
    }
    void Out(const string& name, std::map<int, int>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)ofile << t_.first << " " << t_.second << " ";
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                int a, b;
                ifile >> a >> b;
                t[a] = b;
            }
        }
    }
    void Out(const string& name, std::map<string, string>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.first.size() << " ";
                if (t_.first.size() != 0) ofile << t_.first << " ";
                ofile << t_.second.size() << " ";
                if (t_.second.size() != 0) ofile << t_.second << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                string a(""), b("");
                int nb_;
                ifile >> nb_;
                if (nb_ != 0)ifile >> a;
                ifile >> nb_;
                if (nb_ != 0)ifile >> b;
                t[a] = b;
            }
        }
    }
    void Out(const string& name, std::map<string, int>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.first.size() << " ";
                if (t_.first.size() != 0) ofile << t_.first << " ";
                ofile << t_.second << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                string a("");
                int b, nb_;
                ifile >> nb_;
                if (nb_ != 0)ifile >> a;
                ifile >> b;
                t[a] = b;
            }
        }
    }

    void Out(const string& name, std::vector<std::pair<string, int>>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.first.size() << " ";
                if (t_.first.size() != 0) ofile << t_.first << " ";
                ofile << t_.second << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                string a("");
                int b, nb_;
                ifile >> nb_;
                if (nb_ != 0)ifile >> a;
                ifile >> b;
                t.push_back(std::pair<string, int>(a, b));
            }
        }
    }

    void Out(const string& name, std::map<string, bool>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.first.size() << " ";
                if (t_.first.size() != 0) ofile << t_.first << " ";
                ofile << t_.second << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                string a("");
                bool b;
                int nb_;
                ifile >> nb_;
                if (nb_ != 0)ifile >> a;
                ifile >> b;
                t[a] = b;
            }
        }
    }
    void Out(const string& name, std::map<int, std::tuple<int, int, int, int, int>>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.first << " "
                    << std::get<0>(t_.second) << " "
                    << std::get<1>(t_.second) << " "
                    << std::get<2>(t_.second) << " "
                    << std::get<3>(t_.second) << " "
                    << std::get<4>(t_.second) << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                int a, b0, b1, b2, b3, b4;
                ifile >> a >> b0 >> b1 >> b2 >> b3 >> b4;
                t[a] = std::tuple<int, int, int, int, int>(b0, b1, b2, b3, b4);
            }
        }
    }
    void Out(const string& name, std::map <string, std::map<string, bool>>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << std::endl;
            for (auto& t_ : t)
            {
                string t_str(t_.first);
                Out(name, t_str);
                Out(name, t_.second);
            }
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                string str;
                std::map<string, bool> strb;
                Out(name, str);
                if (str.size() == 0)Functs::MAssert("if (str.size() == 0)");
                Out(name, strb);
                t[str] = strb;
            }
        }
    }

    void Out(const string& name, std::map <string, Vector1i1>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << std::endl;
            for (auto& t_ : t)
            {
                string t_str(t_.first);
                Out(name, t_str);
                Out(name, t_.second);
            }
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                string str;
                Vector1i1 strb;
                Out(name, str);
                if (str.size() == 0)Functs::MAssert("if (str.size() == 0)");
                Out(name, strb);
                t[str] = strb;
            }
        }
    }

    void Out(const string& name, std::map<int, TI3>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.first << " "
                    << std::get<0>(t_.second) << " "
                    << std::get<1>(t_.second) << " "
                    << std::get<2>(t_.second) << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                int a, b0, b1, b2;
                ifile >> a >> b0 >> b1 >> b2;
                t[a] = std::tuple<int, int, int>(b0, b1, b2);
            }
        }
    }
    void Out(const string& name, std::vector<std::map <string, int>>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << std::endl;
            for (auto& t_ : t)
                Out(name, t_);
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = std::vector<std::map <string, int>>(nb, std::map <string, int>());
            for (auto& t_ : t)
                Out(name, t_);
        }
    }
    void Out(const string& name, std::vector<std::map<int, TI3>>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << std::endl;
            for (auto& t_ : t)
                Out(name, t_);
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = std::vector<std::map<int, TI3>>(nb, std::map<int, TI3>());
            for (auto& t_ : t)
                Out(name, t_);
        }
    }
    void Out(const string& name, std::vector<std::tuple<int, int, Vector1i1>>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << std::get<0>(t_) << " " << std::get<1>(t_) << " " << std::get<2>(t_).size()
                    << " " << Functs::Int2String(std::get<2>(t_), false, " ") << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                int a, b, nb_;
                ifile >> a >> b >> nb_;
                Vector1i1 c(nb_, -1);
                for (auto& c_ : c)ifile >> c_;
                t.push_back(std::tuple<int, int, Vector1i1>(a, b, c));
            }
        }
    }
    void Out(const string& name, std::map<std::pair<int, int>, double>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << t_.first.first << " "
                    << t_.first.second << " ";
                ofile << fixed << setprecision(8) << t_.second << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                int a, b; double c;
                ifile >> a >> b >> c;
                t[std::pair<int, int>(a, b)] = c;
            }
        }
    }

    void Out(const string& name, Vector2d& t)
    {
        string tmp, lmp;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " ";
            ofile << fixed << setprecision(8) << t[0] << " " << t[1] << std::endl;
        }
        else
        {
            ifile >> lmp >> tmp >> t[0] >> t[1];
            Read(lmp, tmp, name);
        }
    }


    void Out(const string& name, Vector3d& t)
    {
        string tmp, lmp;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " ";
            ofile << fixed << setprecision(8) << t[0] << " " << t[1] << " " << t[2] << std::endl;;
        }
        else
        {
            ifile >> lmp >> tmp >> t[0] >> t[1] >> t[2];
            Read(lmp, tmp, name);
        }
    }

    void Out(const string& name, glm::dmat4& t)
    {
        string tmp, lmp;

        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " ";
            for (int i = 0; i < 4; i++)for (int j = 0; j < 4; j++) ofile << fixed << setprecision(8) << t[i][j] << " ";
            ofile << std::endl;
        }
        else
        {
            ifile >> lmp >> tmp;
            Read(lmp, tmp, name);
            for (int ii = 0; ii < 4; ii++) for (int jj = 0; jj < 4; jj++)ifile >> t[ii][jj];
        }
    }


    void Out(const string& name, Vector3d2& t)
    {
        string tmp, lmp; int nb;

        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";

            for (auto& t_ : t)
            {
                ofile << t_.size() << " ";
                for (auto& t__ : t_)
                    ofile << fixed << setprecision(8) << t__[0] << " " << t__[1] << " " << t__[2] << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = Vector3d2(nb, Vector3d1());
            for (auto& t_ : t)
            {
                ifile >> nb;
                t_ = Vector3d1(nb, Vector3d());
                for (auto& t__ : t_)
                    ifile >> t__[0] >> t__[1] >> t__[2];
            }
        }
    }


    void Out(const string& name, Vector2d1& t)
    {
        string tmp, lmp; int nb;

        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
                ofile << fixed << setprecision(8) << t_[0] << " " << t_[1] << " ";
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = Vector2d1(nb, Vector2d());
            for (auto& t_ : t)
                ifile >> t_[0] >> t_[1];
        }
    }


    void Out(const string& name, Vector2d2& t)
    {
        string tmp, lmp; int nb;

        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";

            for (auto& t_ : t)
            {
                ofile << t_.size() << " ";
                for (auto& t__ : t_)
                    ofile << fixed << setprecision(8) << t__[0] << " " << t__[1] << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = Vector2d2(nb, Vector2d1());
            for (auto& t_ : t)
            {
                ifile >> nb;
                t_ = Vector2d1(nb, Vector2d());
                for (auto& t__ : t_)
                    ifile >> t__[0] >> t__[1];
            }
        }
    }


    void Out(const string& name, Vector3d1& t)
    {
        string tmp, lmp; int nb;

        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
                ofile << fixed << setprecision(8) << t_[0] << " " << t_[1] << " " << t_[2] << " ";
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            t = Vector3d1(nb, Vector3d());
            for (auto& t_ : t)
                ifile >> t_[0] >> t_[1] >> t_[2];
        }
    }

    void Out(const string& name,
        std::vector<std::pair<Vector3d1, double>>& t)
    {
        if (ol)
        {
            Vector3d2 vecs;
            std::vector<double> ds;
            for (auto& t_ : t)
            {
                vecs.push_back(t_.first);
                ds.push_back(t_.second);
            }
            Out(name, vecs);
            Out(name, ds);
        }
        else
        {
            t.clear();
            Vector3d2 vecs;
            std::vector<double> ds;
            Out(name, vecs);
            Out(name, ds);
            if (vecs.size() != ds.size())Functs::MAssert("vecs.size() != ds.size()");
            for (int i = 0; i < vecs.size(); i++)
            {
                t.push_back(std::pair<Vector3d1, double>(vecs[i], ds[i]));
            }
        }
    }



    void Out(const string& name,
        std::vector<std::pair<Vector3d, Vector3d>>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                ofile << fixed << setprecision(8)
                    << t_.first[0] << " "
                    << t_.first[1] << " "
                    << t_.first[2] << " "
                    << t_.second[0] << " "
                    << t_.second[1] << " "
                    << t_.second[2] << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                Vector3d a, b;
                ifile >> a[0] >> a[1] >> a[2];
                ifile >> b[0] >> b[1] >> b[2];
                t.push_back(std::pair<Vector3d, Vector3d>(a, b));
            }
        }
    }


    void Out(const string& name,
        std::vector<std::tuple<Vector3d, Vector3d, Vector3d, Vector3d>>& t)
    {
        string tmp, lmp; int nb;
        if (ol)
        {
            numb++; ofile << numb << ": "; Check();
            ofile << name << " " << t.size() << " ";
            for (auto& t_ : t)
            {
                //std::get<0>(t_);
                ofile << fixed << setprecision(8)
                    << std::get<0>(t_)[0] << " "
                    << std::get<0>(t_)[1] << " "
                    << std::get<0>(t_)[2] << " "
                    << std::get<1>(t_)[0] << " "
                    << std::get<1>(t_)[1] << " "
                    << std::get<1>(t_)[2] << " "
                    << std::get<2>(t_)[0] << " "
                    << std::get<2>(t_)[1] << " "
                    << std::get<2>(t_)[2] << " "
                    << std::get<3>(t_)[0] << " "
                    << std::get<3>(t_)[1] << " "
                    << std::get<3>(t_)[2] << " ";
            }
            ofile << std::endl;
        }
        else
        {
            t.clear();
            ifile >> lmp >> tmp >> nb;
            Read(lmp, tmp, name);
            for (int i = 0; i < nb; i++)
            {
                Vector3d a, b, c, d;
                ifile >> a[0] >> a[1] >> a[2];
                ifile >> b[0] >> b[1] >> b[2];
                ifile >> c[0] >> c[1] >> c[2];
                ifile >> d[0] >> d[1] >> d[2];
                t.push_back(std::tuple<Vector3d, Vector3d, Vector3d, Vector3d>(a, b, c, d));
            }
        }
    }
};


#endif // _PCE_COMBINE_H

