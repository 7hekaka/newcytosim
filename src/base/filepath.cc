// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "filepath.h"
#include <sys/param.h>
#include <sys/stat.h>
#include <cstdlib>
#include <libgen.h>
#include <unistd.h>
#include <dirent.h>
#include "exceptions.h"


FILE * FilePath::open_file(const char name[], const char mode[2])
{
    if ( name[0] == 0 )
        throw InvalidIO("an empty file name was specified");

    FILE* f = fopen(name, "w");

    if ( !f )
        throw InvalidIO("file `"+std::string(name)+"' could not be opened");
    if ( ferror(f) )
    {
        fclose(f);
        f = nullptr;
        throw InvalidIO("file `"+std::string(name)+"' opened with error");
    }
    return f;
}


std::string FilePath::get_cwd()
{
    char cwd[1024] = { 0 };
    if ( getcwd(cwd, sizeof(cwd)) )
        return cwd;
    else
        return "";
}


bool FilePath::is_file(const char path[])
{
    struct stat s;
    if ( 0 == stat(path, &s) )
        return S_ISREG(s.st_mode);
    return false;
}


bool FilePath::is_dir(const char path[])
{
    struct stat s;
    if ( 0 == stat(path, &s) )
        return S_ISDIR(s.st_mode);
    return false;
}


int FilePath::make_dir(const char name[])
{
    return mkdir(name, S_IRWXU|S_IRWXG|S_IXOTH);
}


int FilePath::change_dir(const char path[])
{
    if ( path && *path )
    {
        if ( chdir(path) )
        {
            perror("Could not change directory");
            errno = 0;
            return -1;
        }
    }
    return 0;
}


int FilePath::change_dir(const char path[], bool make)
{
    int cwd = -1;
    if ( *path )
    {
        cwd = dirfd(opendir("."));
        if ( make )
            (void)mkdir(path, 0777);
        if ( chdir(path) )
        {
            perror("Could not change directory");
            errno = 0;
        }
    }
    return cwd;
}


void FilePath::change_dir(int file)
{
    if ( file >= 0 )
    {
        if ( fchdir(file) )
        {
            perror("Could not change directory");
            errno = 0;
        };
        if ( close(file) )
        {
            perror("Could not close file");
            errno = 0;
        };
    }
}


std::vector<std::string> FilePath::list_dir(const char path[])
{
    std::vector<std::string> res;
    DIR * dp = opendir(path);
    if ( dp )
    {
        struct dirent * ep = readdir(dp);
        while ( ep )
        {
            res.push_back(ep->d_name);
            ep = readdir(dp);
        }
        
        closedir(dp);
    }
    return res;
}


std::vector<std::string> FilePath::list_dir(const char path[], std::string const& ext)
{
    std::vector<std::string> res;
    DIR * dp = opendir(path);
    if ( dp )
    {
        struct dirent * ep = readdir(dp);
        while ( ep )
        {
            std::string name(ep->d_name);
            if ( name.size() > ext.size()
                &&  0==name.compare(name.size()-ext.size(), ext.size(), ext) )
                res.push_back(ep->d_name);
            ep = readdir(dp);
        }
        
        closedir(dp);
    }
    
#if ( 0 )
    for ( size_t i = 0; i < res.size(); ++i )
        std::clog << "   " << res[i] << '\n';
#endif
    return res;
}


std::string FilePath::dir_part(const char path[])
{
    char tmp[PATH_MAX] = { 0 };
    
    if ( !realpath(path, tmp) )
        return ".";
    
    char* res = dirname(tmp);
    
    if ( !res )
        throw InvalidIO("stdlib::dirname() failed");
    
    return std::string(res);
}


std::string FilePath::file_part(const char path[])
{
    char str[MAXPATHLEN] = { 0 };
    strncpy(str, path, MAXPATHLEN);
    char * res = basename(str);
    
    if ( !res )
        throw InvalidIO("stdlib::basename() failed");
    
    return std::string(res);
}


std::string FilePath::full_name(std::string const& dir, std::string const& file)
{
    //if a full path is already specified, we do nothing
    if ( dir.size() > 0  &&  file.size() > 0  &&  file[0] != '/' )
    {
        std::string res = dir + file;
        
        //remove trailling '/' if present
        size_t s = res.size();
        if ( 1 <= s  &&  res[s-1] == '/' )
            res[s-1] = '\0';
        
        return res;
    }
    return file;
}


/*
 This will read the content of the file into buf[]
 size is the amount to which buf was allocated using malloc(), and it can be zero
 If more memory is required, the value of buf[] and size will be updated accordingly
 As `buf` is allocated with 'malloc', it should be eventually released by `free(buf)`
 @returns nullptr in case the file cannot be opened.
 */
void FilePath::read_file(FILE * f, char*& buf, size_t& size)
{
    const size_t chunk = 128;
    if ( !buf || !size )
    {
        size = 8192;
        buf = (char*) malloc(size);
    }
    char * ptr = buf;
    while ( !ferror(f) )
    {
        while ( ptr + chunk >= buf + size )
        {
            ssize_t off = ptr - buf;
            size_t alc = 2 * size;
            buf = (char*) realloc(buf, alc);
            memset(buf+off, 0, alc-off);
            //fprintf(stderr, "%p  %lu ---> %lu\n", &gen, size, alc);
            ptr = buf + off;
            size = alc;
        }
        size_t s = fread(ptr, 1, chunk, f);
        ptr += s;
        if ( s < chunk )
            break;
    }
}


/*
 This will read the content of the file into buf[]
 size is the amount to which buf was allocated using malloc(), and it can be zero
 If more memory is required, the value of buf[] and size will be updated accordingly
 As `buf` is allocated with 'malloc', it should be eventually released by `free(buf)`
 @returns nullptr in case the file cannot be opened.
 */
char* FilePath::read_file(const char filename[], char*& buf, size_t& size)
{
    //fprintf(stderr, "FilePath::read_file(%s)\n", filename);
    FILE * f = fopen(filename, "r");
    if ( f )
    {
        read_file(f, buf, size);
        if ( ferror(f) )
            fprintf(stderr, "Error reading from `%s'\n", filename);
        fclose(f);
        return buf;
    }
    else
    {
        perror(filename);
        errno = 0;
        return nullptr;
    }
}

