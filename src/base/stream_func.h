// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include <iostream>

/// Simple operations on C++ streams
namespace StreamFunc
{
    
    /// remove non-conventional characters
    void clean_stream(std::ostream&, std::istream&);
    
    /// export lines of `val` that are not identical to `ref`
    void diff_stream(std::ostream&, std::istream& val, std::istream& ref);
    
    
    /// copy lines that do not start with character `skip`
    void skip_lines(std::ostream&, std::istream&, char skip, bool verbatim = true);

    /// send lines starting with character `skip` to second stream
    void redirect_lines(std::ostream&, std::ostream&, std::istream&, char sel);

    /// add `prefix` before every line, but skip lines starting with `skip`
    void prefix_lines(std::ostream&, std::istream&, const char prefix[], char verbatim, char skip);

    /// print the line of `istream` indicating the position `pos`, with a marker underline
    void mark_line(std::ostream&, std::istream&, std::streampos pos = -1, const char prefix[]="");

    /// same as `mark_line()`, but output is returned as a string
    std::string marked_line(std::istream&, std::streampos, const char prefix[]);
    
    
    /// extract line containing given `pos`, preserving the current file position
    std::string extract_line(std::istream&, std::streampos pos, size_t& line_nb);
    
    /// extract line containing given `pos`, preserving the current file position
    std::string extract_line(std::istream&, std::streampos pos);
    
    /// extract current line, preserving the current file position
    inline std::string extract_line(std::istream& is) { return extract_line(is, is.tellg()); }

    /// print line located at `pos`, with a line number
    void print_line(std::ostream&, std::istream&, std::streampos pos);
    
    /// extract the lines located between `start` and `end`, with line numbers
    void print_lines(std::ostream&, std::istream&, std::streampos start, std::streampos end, bool all);
    
    /// same as print_lines(), but output is returned as a string
    std::string extract_lines(std::istream&, std::streampos start, std::streampos end);
    
    /// indicate line number and first line
    std::string indicate_lines(std::istream&, std::streampos start, std::streampos end);
    
    /// return number of lines in stream
    size_t number_lines(std::istream&);

    /// return line number corresponding to `pos` (or current position if not specified)
    size_t line_number(std::istream&, std::streampos pos = -1);

    /// replace all occurences of `fnd` by `rep` in `src`. Returns number of replacements done
    std::string replace_string(std::string const& src, std::string const& fnd, std::string const& rep, size_t& cnt);

    /// true if stream has significant unread material
    size_t has_trail(std::istream& is);

}

