#include "libhgp.h"

// libhgp.h includes "cgal.h", which in turn includes "liblgp/liblgp.hpp".
// So liblgp namespace and Functs are available here.

using namespace liblgp;

extern "C" LIBHGP_EXPORT double LGP_Variance_Double(const double* data, int n)
{
    if (!data || n <= 1)
    {
        // For invalid input or too few samples, return 0.0
        return 0.0;
    }

    Vector1d1 vec;
    vec.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
    {
        vec.emplace_back(data[i]);
    }

    return Functs::Variance(vec);
}

extern "C" LIBHGP_EXPORT void LGP_IncreaseVector(int minI, int maxI, Vector1i1 & out)
{
    out = Functs::IncreaseVector(minI, maxI);
}

extern "C" LIBHGP_EXPORT void LGP_ShuffleVector(int size, Vector1i1 & out)
{
    out = Functs::ShuffleVector(size);
}

extern "C" LIBHGP_EXPORT void LGP_Double2String(double value, int p, char* out, int outSize)
{
    if (!out || outSize <= 0)
    {
        return;
    }

    std::string s = Functs::Double2String(value, p);

    // Ensure null-termination and avoid buffer overrun
    if (static_cast<int>(s.size()) >= outSize)
    {
        // Truncate if necessary
        std::string truncated = s.substr(0, static_cast<size_t>(outSize - 1));
        std::memcpy(out, truncated.c_str(), truncated.size());
        out[outSize - 1] = '\0';
    }
    else
    {
        std::memcpy(out, s.c_str(), s.size() + 1); // includes null terminator
    }
}

extern "C" LIBHGP_EXPORT double LGP_String2Double(const char* str)
{
    if (!str)
    {
        return 0.0;
    }
    return Functs::String2Double(std::string(str));
}

extern "C" LIBHGP_EXPORT bool LGP_StringContain(const char* str, const char* sub)
{
    if (!str || !sub)
    {
        return false;
    }
    return Functs::StringContain(std::string(str), std::string(sub));
}

extern "C" LIBHGP_EXPORT void LGP_StringReplace(const char* source,
                                                const char* toReplace,
                                                const char* replaceWith,
                                                char* out,
                                                int outSize)
{
    if (!source || !toReplace || !replaceWith || !out || outSize <= 0)
    {
        return;
    }

    std::string result = Functs::StringReplace(std::string(source),
                                               std::string(toReplace),
                                               std::string(replaceWith));

    if (static_cast<int>(result.size()) >= outSize)
    {
        std::string truncated = result.substr(0, static_cast<size_t>(outSize - 1));
        std::memcpy(out, truncated.c_str(), truncated.size());
        out[outSize - 1] = '\0';
    }
    else
    {
        std::memcpy(out, result.c_str(), result.size() + 1);
    }
}

extern "C" LIBHGP_EXPORT void LGP_WinGetCurDirectory(char* out, int outSize)
{
    if (!out || outSize <= 0)
    {
        return;
    }

    std::string dir = Functs::WinGetCurDirectory();

    if (static_cast<int>(dir.size()) >= outSize)
    {
        std::string truncated = dir.substr(0, static_cast<size_t>(outSize - 1));
        std::memcpy(out, truncated.c_str(), truncated.size());
        out[outSize - 1] = '\0';
    }
    else
    {
        std::memcpy(out, dir.c_str(), dir.size() + 1);
    }
}

extern "C" LIBHGP_EXPORT void LGP_WinGetUserName(char* out, int outSize)
{
    if (!out || outSize <= 0)
    {
        return;
    }

    std::string user = Functs::WinGetUserName();

    if (static_cast<int>(user.size()) >= outSize)
    {
        std::string truncated = user.substr(0, static_cast<size_t>(outSize - 1));
        std::memcpy(out, truncated.c_str(), truncated.size());
        out[outSize - 1] = '\0';
    }
    else
    {
        std::memcpy(out, user.c_str(), user.size() + 1);
    }
}
