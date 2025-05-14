#ifndef GZSTREAM_H
#define GZSTREAM_H

#include <streambuf>
#include <istream>
#include <zlib.h>
#include <stdexcept>

class gzstreambuf : public std::streambuf {
public:
    explicit gzstreambuf(const char* filename, size_t bufsize = 8192)
        : buffer(bufsize), file(gzopen(filename, "rb")) {
        if (!file)
            throw std::runtime_error("Failed to open gzipped file");
        setg(buffer.data(), buffer.data(), buffer.data());
    }

    ~gzstreambuf() {
        if (file) gzclose(file);
    }

protected:
    int_type underflow() override {
        if (gzeof(file)) return traits_type::eof();

        int bytes = gzread(file, buffer.data(), buffer.size());
        if (bytes <= 0) return traits_type::eof();

        setg(buffer.data(), buffer.data(), buffer.data() + bytes);
        return traits_type::to_int_type(*gptr());
    }

private:
    std::vector<char> buffer;
    gzFile file;
};

class gzistream : public std::istream {
public:
    explicit gzistream(const char* filename)
        : std::istream(nullptr), buf(filename) {
        rdbuf(&buf);
    }

private:
    gzstreambuf buf;
};

#endif // GZSTREAM_H
