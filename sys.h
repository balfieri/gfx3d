// Copyright (c) 2017-2022 Robert A. Alfieri
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
// 
// sys.h - system dependencies and utilities
//
// This header provides misc types and functions for:
// - assertions
// - common data types from Model.h
// - strings
// - raw casting between real and uint32_t, or real64 and uint64_t
// - random number generation that is per-thread (and easily implementable in HW)
// - bit twiddling
// - date and time
// - multi-threading 
// - regular expressions
// - networking
//
#ifndef SYSH
#define SYSH

#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <cassert>
#include <vector>
#include <map>
#include <unordered_map>
#include <mutex>
#include <regex>
#include <algorithm>
#include <pthread.h>
#include <time.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netdb.h>

static float uniform( void );  // temporary
#include "Model.h"

// debug
static bool __debug = false;
static bool __rt_debug = false;
static std::mutex __debug_mutex;          // to avoid garble with multiple threads
#define debug_lock() std::lock_guard<std::mutex> guard(__debug_mutex)
#define dout   if (__debug) std::cout
#define rtdout if (__rt_debug || __debug) std::cout
#define wassert(expr) if (!(expr)) std::cout << "WARNING: not true: " << #expr << "\n";
#define dassert(expr, msg) { if ( !(expr) ) die( msg ); }
#define eassert(expr) dassert( expr, "not true: " + std::string(#expr) )

// expose/rename types and constants from Model.h
//
using int128_t = __int128_t;
using uint128_t = __uint128_t;
using real   = Model::real;
using real32 = Model::real32;
using real64 = Model::real64;
constexpr real   REAL_MAX   = MODEL_REAL_MAX;
constexpr real32 REAL32_MAX = MODEL_REAL32_MAX;
constexpr real64 REAL64_MAX = MODEL_REAL64_MAX;
#define REAL_INFINITY    std::numeric_limits<real>::infinity()
#define REAL32_INFINITY  std::numeric_limits<real32>::infinity()
#define REAL64_INFINITY  std::numeric_limits<real64>::infinity()

using real2  = Model::real2;
using vec2   = Model::real2;
using real3  = Model::real3;
using vec3   = Model::real3;
using real3d = Model::real3d;
using vec3d  = Model::real3d;
using real4  = Model::real4;
using vec4   = Model::real4;
using AABB   = Model::AABB;
using aabb   = Model::AABB;
using AABBD  = Model::AABBD;
using aabbd  = Model::AABBD;
using Matrix = Model::Matrix;
using RAY_KIND = Model::RAY_KIND;
class material;
using Ray64  = Model::Ray64;
class ray : public Model::Ray
{
public:
    inline ray() : Model::Ray() {}
    inline ray(const vec3& a, const vec3& b, RAY_KIND kind, real solid_angle=0) : Model::Ray(a, b, kind, solid_angle) {};
    inline ray(const Model::Ray& r) : Model::Ray(r) {};

    inline void copy_stack(const ray& r) 
    {
        top = r.top;
        stack[0] = r.stack[0];
        stack[1] = r.stack[1];
    }

    material* stack[2];
};
using camera = Model::Camera;

// for debugging, set by render_thread() 
static thread_local int my_pixel_x;
static thread_local int my_pixel_y;
static thread_local int my_subpixel_x;
static thread_local int my_subpixel_y;

//--------------------------------------------------------- 
// Strings
//--------------------------------------------------------- 
inline std::string upper( std::string s )
{
    std::string r = "";
    for( uint32_t i = 0; i < s.length(); i++ )
    {
        r += toupper( s[i] );
    }
    return r;
}

//--------------------------------------------------------- 
// Raw Type Casting Between real and uint32_t, or real64 and uint64_t.
//--------------------------------------------------------- 
inline real bits_to_real( uint32_t u ) 
{ 
    real f; 
    memcpy( &f, &u, sizeof(f) ); 
    return f; 
}

inline uint32_t real_to_bits( real r )
{ 
    uint32_t u; 
    memcpy( &u, &r, sizeof(u) ); 
    return u; 
}

inline real64 bits64_to_real64( uint64_t u ) 
{ 
    real64 f; 
    memcpy( &f, &u, sizeof(f) ); 
    return f; 
}

inline uint64_t real64_to_bits64( real64 r )
{ 
    uint64_t u; 
    memcpy( &u, &r, sizeof(u) ); 
    return u; 
}

//--------------------------------------------------------- 
// Random Numbers
//--------------------------------------------------------- 
static thread_local uint32_t my_tid;
static thread_local bool     got_seed = false;
static thread_local uint32_t m_z = 0xbabecafe;       // these are per-thread
static thread_local uint32_t m_w = 0x83417fd1;

inline void register_thread( uint32_t tid )
{
    my_tid = tid; 
    m_z += tid;   // give this thread a unique seed
}

inline void rand_thread_seed( uint64_t seed )
{
    m_z = seed >> 32;
    m_w = seed & 0xffffffffLL;
    got_seed = true;
    //std::cout << "seed=" << seed << "\n";
}

inline uint32_t rand_bits( void )
{
    if ( !got_seed ) {
        std::cout << "ERROR: rand_bits/uniform called before establishing per-thread seed\n";
        exit( 1 );
    }
    m_z = 36969 * (m_z & 65535) + (m_z >> 16);
    m_w = 18000 * (m_w & 65535) + (m_w >> 16);
    uint32_t bits = (m_z << 16) + m_w;
    return bits;
}

inline uint64_t rand_bits64( void )
{
    return (uint64_t( rand_bits() ) << 32) | (uint64_t( rand_bits() ) << 0);
}

inline real uniform( void )
{
    // use frac_w lower bits returned by rand_bits()
    uint32_t rand_frac = rand_bits() & ((1 << MODEL_REAL_FRAC_W)-1);
    double   rand_f    = double(rand_frac) / double(1 << MODEL_REAL_FRAC_W);
    real     rand      = rand_f;
    return rand;
}

inline uint64_t rand_n( uint64_t n )    // returns integer between 0 and n-1
{
    uint64_t bits = rand_bits();
    bits |= uint64_t(rand_bits()) << 32;
    return bits % n;
}

inline bool heads( void )    // returns true or false
{
    return rand_n( 2 ) == 0;
}

//--------------------------------------------------------- 
// Bit Twiddling
//--------------------------------------------------------- 
inline uint32_t bits_reverse( uint32_t bits )
{
    uint32_t r = 0;
    for( uint32_t b = 0; b < 32; b++ )
    {
        uint32_t bit = bits & 1;
        bits >>= 1;
        r |= (bit << (31-b));
    }
    return r;
}

inline uint32_t bits_count_ones( uint32_t bits )
{
    uint32_t r = 0;
    for( uint32_t b = 0; b < 32; b++ )
    {
        if ( bits & 1 ) r++;
        bits >>= 1;
    }
    return r;
}

inline uint32_t bits_count_zeroes( uint32_t bits )
{
    uint32_t r = 0;
    for( uint32_t b = 0; b < 32; b++ )
    {
        if ( !(bits & 1) ) r++;
        bits >>= 1;
    }
    return r;
}

inline void bit_clear( uint32_t& bits, uint32_t b )
{
    bits &= ~(uint32_t(1) << b);
}

inline void bit_set( uint32_t& bits, uint32_t b )
{
    bits |= 1 << b;
}

inline void bit_assign( uint32_t& bits, uint32_t b, bool val )
{
    if ( val ) {
        bit_set( bits, b );
    } else {
        bit_clear( bits, b );
    }
}

inline bool bit_is_zero( uint32_t bits, uint32_t b )
{
    return !((bits >> b) & 1);
}

inline bool bit_is_one( uint32_t bits, uint32_t b )
{
    return (bits >> b) & 1;
}

inline uint32_t bits_insert( uint32_t bits, uint32_t new_bits, uint32_t offset, uint32_t len )
{
    uint32_t mask = ~(0xffffffff << len) << offset;
    bits &= ~mask;
    return bits | (new_bits << offset);
}

inline uint32_t bits_extract( uint32_t bits, uint32_t offset, uint32_t len )
{
    return (bits >> offset) & ((1 << len)-1);
}

inline uint32_t bits_lt( uint32_t b )
{
    b = std::min( b, 32u );
    uint32_t r = 0;
    for( uint32_t i = 0; i < b; i++ )
    {
        r |= 1 << i;
    }
    return r;
}

inline uint32_t bits_le( uint32_t b )
{
    b = std::min( b, 31u );
    uint32_t r = 0;
    for( uint32_t i = 0; i <= b; i++ )
    {
        r |= 1 << i;
    }
    return r;
}

inline uint32_t bits_gt( uint32_t b )
{
    uint32_t r = 0;
    for( uint32_t i = b+1; i < 32; i++ )
    {
        r |= 1 << i;
    }
    return r;
}

inline uint32_t bits_ge( uint32_t b )
{
    uint32_t r = 0;
    for( uint32_t i = b; i < 32; i++ )
    {
        r |= 1 << i;
    }
    return r;
}

inline uint32_t bits_rotate_left( uint32_t bits, uint32_t shift )
{
    shift %= 32;
    return (bits << shift) | (bits >> (32-shift));
}

inline uint32_t bits_rotate_right( uint32_t bits, uint32_t shift )
{
    shift %= 32;
    return (bits >> shift) | (bits << (32-shift));
}

inline uint32_t bits_find_leading_one( uint32_t bits )
{
    for( uint32_t bb = 0; bb < 32; bb++ )
    {
        uint32_t b = 31 - bb;
        if ( (bits >> b) & 1 ) return b;
    }
    die( "bits == 0" );
    return 0;
}

inline uint32_t bits_find_trailing_one( uint32_t bits )
{
    for( uint32_t b = 0; b< 32; b++ )
    {
        if ( (bits >> b) & 1 ) return b;
    }
    die( "bits == 0" );
    return 0;
}

inline int32_t bits_find_leading_one_cbz( uint32_t bits )
{
    return (bits == 0) ? -1 : bits_find_leading_one( bits );
}

inline int32_t bits_find_trailing_one_cbz( uint32_t bits )
{
    return (bits == 0) ? -1 : bits_find_trailing_one( bits );
}

inline uint32_t bits_count_trailing_zeroes( uint32_t bits )
{
    for( uint32_t b = 0; b < 32; b++ )
    {
        if ( (bits >> b) & 1 ) return b;
    }
    return 32;
}

inline uint32_t bits_count_leading_zeroes( uint32_t bits )
{
    for( uint32_t bb = 0; bb < 32; bb++ )
    {
        uint32_t b = 31 - bb;
        if ( (bits >> b) & 1 ) return bb;
    }
    return 32;
}

inline uint32_t bits_find_nth_one_after_with_wrap( uint32_t bits, uint32_t start, uint32_t n )
{
    //--------------------------------------------------------- 
    // We can go around many times until we hit the nth bit.
    // For example, if there's only one '1' and n=10, we'll keep hitting the same bit.
    // 
    // If n == 0, we return 'start' even if bits[start] is a '0'.
    //--------------------------------------------------------- 
    eassert( bits != 0 );
    uint32_t cnt = 0;
    uint32_t i = start;
    for( ; n != 0; n-- )
    {
        do 
        {
            i = (i == 31) ? 0 : (i+1);
        } while( !bit_is_one( bits, i ) ); 
    }
    return i;
}

inline uint32_t bits_find_nth_one_before_with_wrap( uint32_t bits, uint32_t start, uint32_t n )
{
    //--------------------------------------------------------- 
    // Same code as previous, but reverse direction.
    //--------------------------------------------------------- 
    eassert( bits != 0 );
    uint32_t cnt = 0;
    uint32_t i = start;
    for( ; n != 0; n-- )
    {
        do 
        {
            i = (i == 0) ? 31 : (i-1);
        } while( !bit_is_one( bits, i ) ); 
    }
    return i;
}

//--------------------------------------------------------- 
// Date and Time
//
// Times are real64 seconds with high-precision fractional seconds.
// Clock times are since the epoch (1/1/1970).
//--------------------------------------------------------- 
inline real64 clock_time( void ) 
{
    // Return number of real-time seconds since the epoch.
    struct timespec ts;
    clock_gettime( CLOCK_REALTIME, &ts );
    return real64(ts.tv_sec) + real64(ts.tv_nsec)/real64(1000000000);
}

inline void sleep_time( real64 secs ) 
{
    // Sleep for the given seconds which can be fractional.
    // Sleep only if secs is positive.
    // If a signal is caught, this can wake up early, so
    // the caller should always call clock_time() to see
    // what time it really is.
    if ( secs <= 0.0 ) return;
    struct timespec ts;
    ts.tv_sec = secs;
    ts.tv_nsec = (secs - real64(ts.tv_sec)) * 1000000000;
    nanosleep( &ts, NULL );
}

inline void sleep_until_clock_time( real64 until ) 
{
    // calculate time to sleep from current clock time
    real64 secs = until - clock_time();
    sleep_time( secs );
}

//--------------------------------------------------------- 
// Multi-Threading
//--------------------------------------------------------- 
inline uint32_t thread_hardware_core_cnt( void )
{
    return std::thread::hardware_concurrency() / 2;
}

inline uint32_t thread_hardware_thread_cnt( void )
{
    return std::thread::hardware_concurrency();
}

const uint32_t THREAD_CNT_MAX = 100;

using tid_t = pthread_t;

inline void thread_create( tid_t& tid, void *(*fn)(void *), void * arg, size_t stack_size=1024*1024 ) // run with big stack size by default 
{
    pthread_attr_t attr;
    pthread_attr_init( &attr );
    pthread_attr_setstacksize( &attr, stack_size );
    dassert( pthread_create( &tid, &attr, fn, arg ) == 0, "pthread_create() failed" );
}

inline void thread_join( tid_t tid, void** value_ptr=nullptr )
{
    dassert( pthread_join( tid, value_ptr ) == 0, "pthread_join() failed" );
}

struct ThreadInfo {
    tid_t    ptid;
    uint32_t tid;
    uint32_t thread_cnt;
    void    (*fn)(uint32_t, uint32_t, void *);
    void *    arg;
};

void * thread_prestart( void * arg ) 
{
    ThreadInfo * info = reinterpret_cast<ThreadInfo *>( arg );
    info->fn( info->tid, info->thread_cnt, info->arg );
    return nullptr;
}

void thread_parallelize( uint32_t thread_cnt, void (*fn)(uint32_t, uint32_t, void *), void * arg )
{
    if ( thread_cnt == 0 ) thread_cnt = thread_hardware_thread_cnt();
    if ( thread_cnt > THREAD_CNT_MAX ) {
        std::cout << " max number of threads is " << THREAD_CNT_MAX << "\n";
        exit(1);
    }
    //--------------------------------------------------------- 
    // Start other threads (if any).
    // We use Pthreads directly so that we can get a large stack (1 MB)
    // to ensure we never have that kind of memory-trashing problem.
    //--------------------------------------------------------- 
    ThreadInfo * threads = new ThreadInfo[thread_cnt];
    for( uint32_t i = 1; i < thread_cnt; i++ )
    {
        threads[i].tid = i;
        threads[i].thread_cnt = thread_cnt;
        threads[i].fn = fn;
        threads[i].arg = arg;
        thread_create( threads[i].ptid, thread_prestart, &threads[i] );
    }

    //--------------------------------------------------------- 
    // Call the function ourselves with tid=0.
    //--------------------------------------------------------- 
    (*fn)( 0, thread_cnt, arg );

    //--------------------------------------------------------- 
    // Wait for other threads to finish.
    //--------------------------------------------------------- 
    for( uint32_t i = 1; i < thread_cnt; i++ )
    {
        thread_join( threads[i].ptid );
    }

    delete[] threads;
}

//--------------------------------------------------------- 
// Regular Expression Utility Functions
//--------------------------------------------------------- 
inline std::regex regex(std::string re, std::string options="") {
    // validate options
    std::regex::flag_type flags;
    bool got_grammar = false;
    for (size_t i = 0; i < options.length(); i++)
    {
        char ch = options.at(i);
        switch(ch)
        {
            case 'i': flags |= std::regex_constants::icase;                                    break;
            case 'j': flags |= std::regex_constants::ECMAScript;       got_grammar = true;     break;
            case 'p': flags |= std::regex_constants::basic;            got_grammar = true;     break;
            case 'P': flags |= std::regex_constants::extended;         got_grammar = true;     break;
            case 'a': flags |= std::regex_constants::awk;              got_grammar = true;     break;
            case 'g': flags |= std::regex_constants::grep;             got_grammar = true;     break;
            case 'G': flags |= std::regex_constants::egrep;            got_grammar = true;     break;
            default: die( "unknown regex option character: " + std::to_string(ch) );           break;
        }
    }
    if (!got_grammar) flags |= std::regex_constants::ECMAScript;

    // compile regex
    return std::regex(re, flags); 
}

inline bool match(std::string s, const std::regex& regex, std::vector<std::string>& matches) {
    std::smatch sm;
    if (!std::regex_match( s, sm, regex ) ) return false;

    // return matches as a list of strings
    matches = std::vector<std::string>();
    for( size_t i = 0; i < sm.size(); i++ ) matches.push_back(sm[i]);
    return true;
}

inline bool match(std::string s, const std::string re, std::vector<std::string>& matches) {
    return match(s, regex(re), matches);
}

inline bool match(std::string s, const std::string re, std::string options, std::vector<std::string>& matches) {
    return match(s, regex(re, options), matches);
}

inline std::string replace(std::string s, const std::regex& regex, std::string fmt) {
    return std::regex_replace(s, regex, fmt);
}

inline std::string replace(std::string s, std::string re, std::string fmt) {
    return std::regex_replace(s, regex(re), fmt);
}

inline std::string replace(std::string s, std::string re, std::string options, std::string fmt) {
    return std::regex_replace(s, regex(re, options), fmt);
}

inline std::string indent_str(uint32_t cnt) {
    std::string s = "";
    for (uint32_t i=0; i < cnt; i++) s += " ";
    return s;
}

//--------------------------------------------------------- 
// Networking Utility Functions
//
// Work with both IPv4 or IPv6.
//--------------------------------------------------------- 

inline std::string errno_str( void )
{
    char s[256];
    strerror_r( errno, s, sizeof(s) );
    return std::to_string( errno ) + " (" + s + ")";
}

using socket_id_t = int;
using socket_addr_info_t = struct addrinfo;
using socket_addr_t = struct sockaddr;
using socket_addrlen_t = socklen_t;

socket_addr_info_t * socket_addr_info_alloc( std::string ip_addr, uint32_t port )
{
    struct addrinfo hints;
    memset( &hints, 0, sizeof(hints) );
    hints.ai_family = AF_UNSPEC; // can use IPv4 or IPv6
    hints.ai_socktype = SOCK_DGRAM;
    std::string service = std::to_string( port );
    socket_addr_info_t * addr_info_ptr;
    int ret = getaddrinfo( ip_addr.c_str(), service.c_str(), &hints, &addr_info_ptr );
    dassert( ret == 0, "getaddrinfo() failed for socket_addr_create()" );
    return addr_info_ptr;
}

void socket_addr_info_free( socket_addr_info_t * info )
{
    freeaddrinfo( info );
}

// ip_addr == "" ===> broadcast
socket_addr_info_t * socket_addr_info_udp_alloc( std::string ip_addr, uint32_t port, bool is_ipv4=true )
{
    uint32_t family = is_ipv4 ? AF_INET : AF_INET6;
    socket_addr_t * addr;
    socket_addrlen_t addr_len;
    if ( family == AF_INET ) {
        struct sockaddr_in * addr_in = new struct sockaddr_in;
        addr                     = reinterpret_cast<struct sockaddr *>( addr_in );
        addr_len                 = sizeof( struct sockaddr_in );
        addr_in->sin_family      = family;
        addr_in->sin_port        = htons( port );
        if ( ip_addr == "" ) {
            addr_in->sin_addr.s_addr = INADDR_ANY;
        } else {
            inet_pton( family, ip_addr.c_str(), &addr_in->sin_addr );
        }
        memset( addr_in->sin_zero, 0, sizeof(addr_in->sin_zero) );
    } else {
        struct sockaddr_in6 * addr_in = new struct sockaddr_in6;
        addr                     = reinterpret_cast<struct sockaddr *>( addr_in );
        addr_len                 = sizeof( struct sockaddr_in6 );
        addr_in->sin6_len        = addr_len;
        addr_in->sin6_family     = AF_INET6;
        addr_in->sin6_flowinfo   = 0;
        addr_in->sin6_family     = family;
        addr_in->sin6_port       = htons( port );
        if ( ip_addr == "" ) {
            addr_in->sin6_addr   = in6addr_any;
        } else {
            inet_pton( family, ip_addr.c_str(), &addr_in->sin6_addr );
        }
    }

    char canonname[1] = "";
    socket_addr_info_t * info = new socket_addr_info_t;
    info->ai_flags     = 0;
    info->ai_family    = family;
    info->ai_socktype  = SOCK_DGRAM;
    info->ai_protocol  = 0;
    info->ai_addr      = addr;
    info->ai_addrlen   = addr_len;
    info->ai_canonname = canonname;
    info->ai_next      = nullptr;

    return info;
}

void socket_addr_info_udp_free( socket_addr_info_t * info )
{
    if ( info->ai_family == AF_INET ) {
        struct sockaddr_in * addr_in = reinterpret_cast<struct sockaddr_in *>( info->ai_addr );
        delete addr_in;
    } else {
        dassert( info->ai_family == AF_INET6, "ai_family is not AF_INET or AF_INET6" );
        struct sockaddr_in6 * addr_in = reinterpret_cast<struct sockaddr_in6 *>( info->ai_addr );
        delete addr_in;
    }
    info->ai_addr = nullptr;
    delete info;
}

void udp_socket_create( socket_id_t& sid, socket_addr_t& local_addr, socket_addrlen_t& local_addr_len, const socket_addr_info_t * local_addr_info, bool non_blocking=true )
{
    // try each possible local addr in the list
    for( const socket_addr_info_t * info = local_addr_info; info != nullptr; info = info->ai_next )
    {
        sid = socket( info->ai_family, info->ai_socktype, info->ai_protocol );
        if ( sid < 0 ) continue;

        int ret;
        bool is_broadcast;
        if ( info->ai_family == AF_INET ) {
            struct sockaddr_in  * addr_in  = reinterpret_cast<struct sockaddr_in  *>( info->ai_addr );
            is_broadcast = addr_in->sin_addr.s_addr == INADDR_ANY;
        } else {
            dassert( info->ai_family == AF_INET6, "ai_family is not AF_INET or AF_INET6" );
            struct sockaddr_in6 * addr_in6 = reinterpret_cast<struct sockaddr_in6 *>( info->ai_addr );
            is_broadcast = true;
            for( uint32_t i = 0; i < 16; i++ )
            {
                is_broadcast &= addr_in6->sin6_addr.s6_addr[i] == in6addr_any.s6_addr[i];
            }
        }
        if ( is_broadcast ) {
            // This is required for broadcast setup.
            // See: https://www.cs.ubbcluj.ro/~dadi/compnet/labs/lab3/udp-broadcast.html
            int so_broadcast = 1;
            ret = setsockopt( sid, SOL_SOCKET, SO_BROADCAST, &so_broadcast, sizeof(so_broadcast) );
            dassert( ret == 0, "setsockopt() failed for udp_socket_create() with broadcast addr errno=" + errno_str() );
        }

        ret = bind( sid, info->ai_addr, info->ai_addrlen );
        if ( ret < 0 ) {
            close( sid );
            continue;
        }

        local_addr     = *local_addr_info->ai_addr; 
        local_addr_len = local_addr_info->ai_addrlen;
        if ( non_blocking ) {
            ret = fcntl( sid, F_SETFL, fcntl( sid, F_GETFL ) | O_NONBLOCK );
            dassert( ret == 0, "fcntl() failed for udp_socket_create() errno=" + errno_str() );
        }
        return; 
    }
    die( "udp_socket_create() could not find any bindable socketaddr errno=" + errno_str() );
}

void udp_socket_create_unicast( socket_id_t& sid, socket_addr_t& local_addr, socket_addrlen_t& local_addr_len, std::string ip_addr, uint32_t port, bool non_blocking=true, bool is_ipv4=true )
{
    socket_addr_info_t * info = socket_addr_info_udp_alloc( ip_addr, port, is_ipv4 );
    udp_socket_create( sid, local_addr, local_addr_len, info, non_blocking );
    socket_addr_info_udp_free( info );
}

void udp_socket_create_broadcast( socket_id_t& sid, socket_addr_t& local_addr, socket_addrlen_t& local_addr_len, uint32_t port, bool non_blocking=true, bool is_ipv4=true )
{
    socket_addr_info_t * info = socket_addr_info_udp_alloc( "", port, is_ipv4 );
    udp_socket_create( sid, local_addr, local_addr_len, info, non_blocking );
    socket_addr_info_udp_free( info );
}

void udp_socket_destroy( socket_id_t sid )
{
    close( sid );
}

void udp_socket_recvfrom( size_t& byte_cnt, socket_id_t sid, void * buffer, size_t buffer_len, socket_addr_t& remote_addr, socket_addrlen_t& remote_addr_len, std::string * remote_ip_addr_ptr=nullptr )
{
    remote_addr_len = sizeof( socket_addr_t );
    int ret = recvfrom( sid, buffer, buffer_len, 0, &remote_addr, &remote_addr_len );
    if ( ret < 0 ) {
        dassert( errno == EAGAIN || errno == EWOULDBLOCK, "recvfrom() failed for udp_socket_recvfrom() errno=" + errno_str() );
        byte_cnt = 0;
    } else {
        byte_cnt = ret;
        if ( remote_ip_addr_ptr != nullptr ) {
            // return string form of remote IPv4 or IPv6 address
            void * ip_addr;
            if ( remote_addr.sa_family == AF_INET ) {
                struct sockaddr_in * addr_in = reinterpret_cast<struct sockaddr_in *>( &remote_addr );
                ip_addr = &addr_in->sin_addr;
            } else {
                dassert( remote_addr.sa_family == AF_INET6, "remote_addr.sa_family is not AF_INET or AF_INET6, got " + std::to_string(remote_addr.sa_family) );
                struct sockaddr_in6 * addr_in6 = reinterpret_cast<struct sockaddr_in6 *>( &remote_addr );
                ip_addr = &addr_in6->sin6_addr;
            }
            char ip_addr_cs[INET6_ADDRSTRLEN];
            const char * ret = inet_ntop( remote_addr.sa_family, ip_addr, ip_addr_cs, sizeof(ip_addr_cs) );
            dassert( ret != nullptr, "inet_ntop() failed errno=" << errno_str() );
            *remote_ip_addr_ptr = ip_addr_cs;
        }
    }
}

void udp_socket_sendto( size_t& byte_cnt, socket_id_t sid, void * buffer, size_t buffer_len, const socket_addr_t& local_addr, size_t local_addr_len )
{
    int ret = sendto( sid, buffer, buffer_len, 0, &local_addr, local_addr_len );
    if ( ret < 0 ) {
        dassert( errno == EAGAIN || errno == EWOULDBLOCK, "sendto() failed for udp_socket_sendto() errno=" + errno_str() );
        byte_cnt = 0;
    } else {
        byte_cnt = ret;
    }
}

#endif
