#ifndef INVALID_OPERATION_EXCEPTION_H__
#define INVALID_OPERATION_EXCEPTION_H__

#pragma once

#include <stdexcept>

class InvalidOperationException : public std::exception
{
private:
    std::string msg;

public:
    InvalidOperationException(const std::string& message = "") : msg(message)
    {
    }

    const char * what() const throw()
    {
        return msg.c_str();
    }
};

#endif
