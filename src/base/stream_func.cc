// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "stream_func.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cctype>
#include <algorithm>

void StreamFunc::clean_stream(std::ostream& os, std::istream& is)
{
    char c = 0;
    while ( is.good() )
    {
        is.get(c);
        
        // terminate the line for new-line and cariage-return
        if ( c == '\r' )
            os.put('\n');
        // all type of spaces are substituted
        else if ( isspace(c) )
            os.put(' ');
        // non=printable characters are removed
        else if ( isprint(c) )
            os.put(c);
        else
            std::cerr << "unprintable ascii "<< (int)c << " found\n";
    }
}


void StreamFunc::diff_stream(std::ostream& os, std::istream& val, std::istream& ref)
{
    std::string val_l, ref_l;
    val.seekg(0);
    ref.seekg(0);
    
    while ( val.good() )
    {
        std::getline(val, val_l);
        std::getline(ref, ref_l);
#if ( 0 )
        // print any line containing '{' or '}' 
        bool par = ( std::string::npos != ref_l.find_first_of("{}") );
        if ( val_l != ref_l || par )
#else
        if ( val_l != ref_l )
#endif
        {
            os << val_l << '\n';
            //os << val_l << "  (" << ref_l << ")\n";
        }
    }
}


void StreamFunc::skip_lines(std::ostream& os, std::istream& is, char skip)
{
    std::string line;
    while ( is.good() )
    {
        std::getline(is, line);
        if ( is.fail() )
            break;
        if ( line.empty() || line[0] != skip )
            os << line << '\n';
    }
}


void StreamFunc::redirect_lines(std::ostream& os, std::ostream& alt, std::istream& is, char sel)
{
    std::string line;
    
    while ( is.good() )
    {
        std::getline(is, line);
        if ( is.fail() )
            break;
        if ( line.size() > 0 && line[0] == sel )
        {
            alt << line;
            if ( !is.eof() )
                alt << '\n';
        }
        else
        {
            os << line;
            if ( !is.eof() )
                os << '\n';
        }
    }
}


void StreamFunc::prefix_lines(std::ostream& os, std::istream& is, const char prefix[],
                              char verbatim, char skip)
{
    std::string line;
    
    while ( is.good() )
    {
        std::getline(is, line);
        if ( is.fail() )
            break;
        if ( line.size() < 1 )
            os << '\n';
        else if ( line[0] == skip )
            ;
        else if ( line[0] == verbatim )
        {
            os << line;
            if ( !is.eof() )
                os << '\n';
        }
        else
        {
            os << prefix << line;
            if ( !is.eof() )
                os << '\n';
        }
    }
}


/**
 The alignment of the vertical bar should match the one in PREF
 */
static void print_first_line(std::ostream& os, size_t num, std::string const& line)
{
    os << " in" << std::setw(6) << num << " | " << line << '\n';
}

/**
 The alignment of the vertical bar should match the one in PREF
 */
static void print_other_line(std::ostream& os, size_t num, std::string const& line)
{
    os << std::setw(9) << num << " | " << line << '\n';
}

/**
 The alignment of the vertical bar should match the one in PREF
 */
static void print_other_line(std::ostream& os, const char prefix[], std::string const& line)
{
    if ( prefix && *prefix )
        os << prefix << " " << line << '\n';
    else
        os << line << '\n';
}


/**
 Print the one line extracted from `is` containing `pos` and indicate
 the position `pos` with a arrowhead in a second line below.
 */
void StreamFunc::mark_line(std::ostream& os, std::istream& is, std::streampos pos, const char prefix[])
{
    is.clear();
    std::streampos sos = is.tellg(), isp = sos;
    is.seekg(0);

    // get the line containing 'pos'
    std::string line;
    while ( is.good()  &&  is.tellg() <= pos )
    {
        isp = is.tellg();
        std::getline(is, line);
    }
    std::streamoff off = pos - isp;
    // reset stream
    is.clear();
    is.seekg(sos);

    std::string sub;
    std::streamoff i = 0;
    while ( i < off )
    {
        int c = line[i++];
        if ( c == 0 )
            break;
        sub.push_back(isspace(c)?(char)c:' ');
    }
    sub.push_back('^');
    //sub.append(" ("+std::string(1, is.peek())+")");
    os << '\n' << prefix << " " << line;
    os << '\n' << prefix << " " << sub;
}


void StreamFunc::mark_line(std::ostream& os, std::istream& is)
{
    is.clear();
    mark_line(os, is, is.tellg(), nullptr);
}


std::string StreamFunc::marked_line(std::istream& is, std::streampos pos, const char prefix[])
{
    std::ostringstream oss;
    mark_line(oss, is, pos, prefix);
    return oss.str();
}


/**
 Output enough lines to cover the area specified by [start, end].
 Each line is printed in full and preceded with a line number
 */
void StreamFunc::print_lines(std::ostream& os, std::istream& is,
                             std::streampos start, std::streampos end)
{
    if ( !is.good() )
        is.clear();
    
    std::streampos isp = is.tellg();
    is.seekg(0);
    std::string line;
    
    size_t cnt = 0;
    while ( is.good()  &&  is.tellg() <= start  )
    {
        std::getline(is, line);
        ++cnt;
    }

    print_first_line(os, cnt, line);
    while ( is.good() &&  is.tellg() < end )
    {
        std::getline(is, line);
        ++cnt;
        if ( !std::all_of(line.begin(),line.end(),isspace) )
            print_other_line(os, cnt, line);
    }

    is.clear();
    is.seekg(isp);
}


std::string StreamFunc::extract_lines(std::istream& is, std::streampos s, std::streampos e)
{
    std::ostringstream oss;
    print_lines(oss, is, s, e);
    return oss.str();
}


std::string StreamFunc::extract_line(std::istream& is, std::streampos pos, size_t& cnt)
{
    if ( !is.good() )
        is.clear();
    std::streampos isp = is.tellg();

    is.seekg(0);
    std::string line;
    cnt = 0;
    
    while ( is.good()  &&  is.tellg() <= pos )
    {
        std::getline(is, line);
        ++cnt;
    }
    
    is.clear();
    is.seekg(isp);
    return line;
}


std::string StreamFunc::extract_line(std::istream& is, std::streampos pos)
{
    size_t cnt;
    return extract_line(is, pos, cnt);
}


void StreamFunc::print_line(std::ostream& os, std::istream& is, std::streampos pos)
{
    size_t cnt;
    std::string line = extract_line(is, pos, cnt);
    os << "line" << std::setw(6) << cnt << " | " << line << '\n';
}


size_t StreamFunc::number_lines(std::istream& is)
{
    if ( !is.good() )
        is.clear();
    
    std::streampos isp = is.tellg();
    is.seekg(0);
    
    size_t cnt = 0;
    std::string line;
    
    while ( is.good() )
    {
        std::getline(is, line);
        ++cnt;
    }
    
    is.clear();
    is.seekg(isp);
    return cnt;
}


size_t StreamFunc::line_number(std::istream& is, std::streampos pos)
{
    if ( !is.good() )
        is.clear();
    
    std::streampos isp = is.tellg();
    is.seekg(0);
    
    if ( pos == -1 )
        pos = isp;
    
    size_t cnt = 0;
    std::string line;
    
    while ( is.good()  &&  is.tellg() <= pos )
    {
        std::getline(is, line);
        ++cnt;
    }
    
    is.clear();
    is.seekg(isp);
    return cnt;
}


std::string StreamFunc::replace_string(std::string const& src, std::string const& fnd,
                                       std::string const& sub, size_t& cnt)
{
    cnt = 0;
    std::string res;
    std::string::size_type P = 0;
    std::string::size_type Q = src.find(fnd, 0);
    res.append(src, 0, Q);
    while ( Q != std::string::npos )
    {
        ++cnt;
        res.append(sub);
        P = Q + fnd.size();
        Q = src.find(fnd, P);
        res.append(src, P, Q-P);
    }
    return res;
}


/// return `true` if stream contains unread non-space character(s)
bool StreamFunc::has_trail(std::istream& is)
{
    int c = is.get();
    while ( isspace(c) )
        c = is.get();
    if ( c != EOF )
    {
        is.unget();
        return true;
    }
    return false;
}

