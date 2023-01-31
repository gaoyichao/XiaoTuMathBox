#ifndef XTMB_EXCEPTION_H
#define XTMB_EXCEPTION_H

#include <stdint.h>
#include <cassert>
#include <iostream>


namespace xiaotu {
namespace math {

    class NotDeepCopyException : public std::exception {
        public:
            NotDeepCopyException(std::string const & msg)
                : _msg(msg)
            {

            }

            NotDeepCopyException(char const * msg)
                : _msg(msg)
            {

            }

            const char * what() const noexcept
            {
                return _msg.c_str();
            }

        private:
            std::string _msg;
    };


}
}

#endif
