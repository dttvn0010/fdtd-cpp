#ifndef __BASE_H__
#define __BASE_H__
#include <cstdlib>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <algorithm>

#include <type_traits>
#include <typeinfo>
#include <cstdarg>
#include <iostream>
#include <initializer_list>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <tuple>
#include <thread>
#include <memory>

#define isize int64_t
#define DEBUG
#define Using(...) if(__VA_ARGS__; true)
#define RBlock(statements) [&] statements ()

template<typename T>
concept Clonable = std::move_constructible<T> && requires(const T t)
{
    {t.clone()} -> std::convertible_to<T>;
};

template<typename T>
concept Printable = requires(std::ostream & os, const T t)
{
    os << t;
};

inline static void panic(const char* fmt, ...) {
    va_list arglist;

    va_start(arglist, fmt);
    int len = vsnprintf(NULL, 0, fmt, arglist);
    va_end(arglist);

    char* tmp = (char*)malloc(len + 1);

    va_start(arglist, fmt);
    vsnprintf(tmp, len + 1, fmt, arglist);
    va_end(arglist);

    std::string msg = tmp;
    free(tmp);
    std::cout << "Error: " + msg << std::endl;
    throw msg;
}


template<typename T>
struct Optional {
    T value = T();
    bool is_null = true;

    Optional() {}

    Optional(const T& value)
    {
        this->value = value;
        this->is_null = false;
    }

    Optional(T&& value)
    {
        this->value = std::move(value);
        this->is_null = false;
    }

    static Optional none()
    {
        return Optional();
    }

};
typedef Optional<int> Integer;

class Object
{
protected:
    size_t _signature = rand() + 1;

public:
    size_t get_signature() const
    {
        return _signature;
    }
    ~Object()
    {
        _signature = 0;
    }
};

template<typename T>
concept IsObject = std::derived_from<T, Object>;

template<typename  T>
class Ptr
{
    size_t _signature = 0;
    T* _data = NULL;

public:
    Ptr()
    {
    }

    Ptr(T* _data)
    {
        this->_data = _data;
        if(_data)
        {
            _signature = ((Object*)_data)->get_signature();
        }
        else
        {
            _signature = 0;
        }
    }

    T* get_unsafe()
    {
        return _data;
    }

    const T* get_unsafe() const
    {
        return _data;
    }

    T* get()
    {

        if (!_data || _signature != ((Object*)_data)->get_signature())
        {
            throw std::string("Object is already destroyed\n");
        }

        return _data;
    }

    const T* get() const
    {

        if (!_data || _signature != ((Object*)_data)->get_signature())
        {
            throw std::string("Object is already destroyed\n");
        }

        return _data;
    }

    T* operator -> ()
    {
        return get();
    }

    const T* operator -> () const
    {
        return get();
    }

    T& operator*()
    {
        return *(operator->());
    }

    const T& operator*() const
    {
        return *(operator->());
    }
};

template <const char end = '\n'>
void print() {
    std::cout << end;
}

template <const char end = '\n', typename T>
void print(const T& t)
{
    std::cout << t << end;
}

template<const char end = '\n', const char sep = ' ', typename T, typename... Args>
void print(const T& t, const Args&... args) // recursive variadic function
{
    std::cout << t << sep;
    print(args...);
}

template<typename T>
class List;

template<typename ITR>
class IteratorRange
{
public:
    ITR begin;
    ITR end;

    IteratorRange(ITR begin, ITR end)
    {
        this->begin = begin;
        this->end = end;
    }
};

template<typename ITR>
class ParallelIterator
{
    std::vector<IteratorRange<ITR>> ranges;

public:
    ParallelIterator(const std::vector<IteratorRange<ITR>>& ranges)
    {
        this->ranges = ranges;
    }

    ParallelIterator(std::vector<IteratorRange<ITR>>&& ranges)
    {
        this->ranges = ranges;
    }

    template<typename Fn>
    void for_each(const Fn& fn)
    {
        std::vector<std::thread> threads;
        for (auto& r : ranges)
        {
            threads.push_back(std::thread(
                [&fn](IteratorRange<ITR>* r) {
                    auto ptr = r->begin;
                    while (ptr != r->end)
                    {
                        fn(*ptr);
                        ++ptr;
                    }
                }
            , &r));
        }

        for (auto& thr : threads)
        {
            thr.join();
        }
    }

    template<typename Fn>
    auto map(const Fn& fn)
    {
        using U = decltype(fn(*ranges[0].begin));
        int n_thread = ranges.size();
        List<List<U>> results(n_thread);

        auto chunk_fn = [&fn](IteratorRange<ITR>* r, List<U>* result) {
            auto ptr = r->begin;
            while (ptr != r->end)
            {
                result->append(fn(*ptr));
                ++ptr;
            }
        };

        std::vector<std::thread> threads;
        for (int i = 0; i < n_thread; i++)
        {
            threads.push_back(std::thread(chunk_fn, &ranges[i], &results[i]));
        }

        for (auto& thr : threads)
        {
            thr.join();
        }

        List<U> lst;
        for (auto& result : results)
        {
            lst.insert(lst.end(), result.begin(), result.end());    // concatenate
        }
        return lst;
    }
};

template<typename COLLECTION>
auto get_par_iter(COLLECTION& lst, int n_threads)
{
    using ITR = decltype(lst.begin());

    std::vector<IteratorRange<ITR>> iterators;
    double start = 0;
    double step = ((double)lst.size()) / n_threads;

    for (int i = 0; i < n_threads; i++)
    {
        double stop = start + step;
        auto begin = lst.begin();
        begin += (int)(start + 0.5);
        auto end = lst.begin();
        end += (int)(stop + 0.5);
        iterators.push_back(IteratorRange(begin, end));
        start = stop;
    }
    return ParallelIterator(iterators);
}

template<typename T>
class ListView;

template<typename T>
class ListItem
{
    T* _data = NULL;
    List<T>* _parent = NULL;
public:
    ListItem()
    {
    }

    ListItem(List<T>& lst, int idx)
    {
        _parent = &lst;
        _data = &lst.at(idx);
        _parent->refcnt += 1;
    }

    ListItem(const List<T>& lst, int idx)
    {
        _parent = &lst;
        _data = &lst.at(idx);
        _parent->refcnt += 1;
    }

    ~ListItem()
    {
        if(_parent) _parent->refcnt -= 1;
    }

    T* get() const
    {
        return _data;
    }

    T* operator -> () const
    {
        return get();
    }

    T& operator*() const
    {
        return *(operator->());
    }
};

template<typename T>
class List : public Object, public std::vector<T>
{
    T* get_data()
    {
        return std::vector<T>::data();
    }

    const T* get_data() const
    {
        return std::vector<T>::data();
    }

    List(const List& lst) = delete;
    void operator=(const List&) = delete;
    template<typename K, typename V> friend class Dict;
    friend class ListItem<T>;
    int refcnt = 0;

    void pushback(const T& item)
    {
        std::vector<T>::push_back(item);
    } 

    void pushback(T&& item)
    {
        std::vector<T>::push_back(item);
    } 

public:
    List() : std::vector<T>() {}

    List(const std::vector<T>& v) : std::vector<T>(v) {}

    List(int sz) : std::vector<T>(sz) {}

    List(std::vector<T>&& v) : std::vector<T>(v) {}

    List(List&& lst) : std::vector<T>((std::vector<T>&&)lst) {}

    List(const std::initializer_list<T>& items) requires std::copy_constructible<T> : 
        std::vector<T>(items){}

    List(const std::initializer_list<T>& items) requires Clonable<T>
        : std::vector<T>()
    {
        for (auto& item : items) {
            std::vector<T>::push_back(item.clone());
        }
    }

    List clone() const requires std::copy_constructible<T> 
    {
        return List((const std::vector<T>&) (*this));
    }

    List clone() const requires Clonable<T> 
    {
        List<T> lst;
        for (auto& item : *this) {
            lst.push_back(item.clone());
        }
        return lst;
    }

    void operator=(List&& lst)
    {
        std::vector<T>::operator=((std::vector<T>&&)lst);
    }

    void append(const T& val)
    {
        this->push_back(val);
    }

    void append(T&& val)
    {
        this->push_back(std::move(val));
    }

    ListItem<T> item_at(int idx)
    {
        return ListItem(*this, idx);
    }

    ListItem<const T> item_at(int idx) const
    {
        return ListItem(*this, idx);
    }

    ListItem<T> first()
    {
        if(size() == 0) return ListItem<T>();
        return ListItem<T>(*this, 0);
    }

    ListItem<const T> first() const
    {
        if(size() == 0) return ListItem<T>();
        return ListItem<T>(*this, 0);
    }

    ListItem<T> last()
    {
        if(size() == 0) return ListItem<T>();
        return ListItem<T>(*this, (int) size() - 1);
    }

    ListItem<const T> last() const
    {
        if(size() == 0) return ListItem<T>();
        return ListItem<T>(*this, (int) size() - 1);
    }

    ListView<const T> slice(size_t start, size_t stop, int step = 1) const
    {
        assert(start >= 0 && start <= stop && stop <= size() && step > 0);
        return ListView<const T>(
            get_data() + start,
            (int)stop - start,
            step,
            &_signature
        );
    }

    ListView<T> slice(size_t start, size_t stop, int step = 1)
    {
        assert(start >= 0 && start <= stop && stop <= size() && step > 0);
        return ListView<T>(
            get_data() + start,
            (int)stop - start,
            step,
            &_signature
        );
    }

    ListView<T> as_view()
    {
        return ListView(get_data(), std::vector<T>::size(), 1, &_signature);
    }

    ListView<const T> as_view() const
    {
        return ListView(get_data(), std::vector<T>::size(), 1, &_signature);
    }

    size_t size() const
    {
        return std::vector<T>::size();
    }

    template<typename Func>
    void sort(Func func)
    {
        std::sort(std::vector<T>::begin(), std::vector<T>::end(), func);
    }

    auto par_iter(int n_threads)
    {
        return get_par_iter(*this, n_threads);
    }

    auto par_iter(int n_threads) const
    {
        return get_par_iter(*this, n_threads);
    }
};


template<typename T>
class ListViewIterator {
    int _step;
    T* _data;
public:

    ListViewIterator(T* data, int step = 1) {
        _data = data;
        _step = step;
    }

    ListViewIterator& operator++() {
        _data += _step;
        return *this;
    }

    bool operator != (ListViewIterator& rhs) {
        return _data != rhs._data;
    }

    T& operator*() const {
        return *_data;
    }
};

template<typename T>
class ListView
{
    int _size;
    int _step;
    T* _data;
    const size_t* _p_signature;
    friend class List<T>;
    friend class List<typename std::remove_cv<T>::type>;

    ListView(T* data, int size, int step, const size_t* p_signature)
    {
        _data = data;
        _size = size;
        _step = step;
        _p_signature = p_signature;
    }

    void check_valid()
    {
#ifdef DEBUG
        if (*_p_signature != typeid(List<T>).hash_code())
        {
            printf("List is invalid\n");
        }
#endif
    }
public:
    const T& operator[](int idx_) const {
        check_valid();
        int idx = idx < 0 ? _size + idx_ : idx_;
#ifdef DEBUG
        if (idx < 0 || idx >= _size)
        {
            panic("Index %d is out of range %d\n", idx_, _size);
        }
#endif

        return _data[idx];
    }

    T& operator[](int idx_) {
        check_valid();
        int idx = idx < 0 ? _size + idx_ : idx_;
#ifdef DEBUG
        if (idx < 0 || idx >= _size)
        {
            panic("Index %d is out of range %d\n", idx_, _size);
        }
#endif

        return _data[idx];
    }

    ListView slice(int start, int stop, int step = 1)
    {
        assert(start >= 0 && start <= stop && stop <= size() && step > 0);
        return ListView(
            _data + start * _step,
            stop - start,
            step * _step,
            _p_signature
        );
    }

    ListViewIterator<T> begin() const
    {
        return ListViewIterator(_data, _step);
    }

    ListViewIterator<T> end() const
    {
        return ListViewIterator(_data + _step * _size, _step);
    }

    int size() const
    {
        return _size;
    }

    List<T> as_list() const
    {
        List<T> lst(_size);
        for (int i = 0; i < _size; i++)
        {
            lst[i] = _data[i];
        }
        return lst;
    }
};

template<typename  T>
class ListObject : public List<std::unique_ptr<T>>
{
public:
    void append(T&& item)
    {
        List<std::unique_ptr<T>>::append(std::make_unique<T>(std::move(item)));
    }

    size_t size() const 
    {
        return this->size();
    }

    Ptr<T> operator[](int idx)
    {
        return this->at(idx).get();
    }

    Ptr<const T> operator[](int idx) const
    {
        return this->at(idx).get();
    }
};

template<IsObject T>
class ListSharedObject : public List<std::shared_ptr<T>>
{
public:
    void append(T&& item)
    {
        List<std::unique_ptr<T>>::append(std::make_shared<T>(std::move(item)));
    }

    void append(const T& item)
    {
        List<std::unique_ptr<T>>::append(std::make_shared<T>(item));
    }

    size_t size() const 
    {
        return this->size();
    }

    Ptr<T> operator[](int idx)
    {
        return this->at(idx).get();
    }

    Ptr<const T> operator[](int idx) const
    {
        return this->at(idx).get();
    }
};

template<IsObject T>
using ListPtr = List<Ptr<T>>;

template<typename LIST>
static std::ostream& print_list(std::ostream& os, const LIST& lst) {
    os << "[";
    int i = 0, sz = lst.size();
    for (const auto& it : lst) {
        os << it;
        if (++i < sz) os << ", ";
    }
    os << "]";
    return os;
}

template<Printable T>
static std::ostream& operator<<(std::ostream& os, const List<T>& lst)
{
    return print_list(os, lst);
}

template<Printable T>
static std::ostream& operator<<(std::ostream& os, const ListView<T>& lst)
{
    return print_list(os, lst);
}

template<typename T>
class Set final : public Object, public std::unordered_set<T>
{
    Set(const Set& s) = delete;
    void operator=(const Set&) = delete;
    friend class List<Set<T>>;
    template<typename K, typename V> friend class Dict;
    template<typename K, typename V> friend struct std::pair;
public:
    Set() : std::unordered_set<T>() {
        
    }

    Set(const std::unordered_set<T>& s) : std::unordered_set<T>(s) {}

    Set(std::unordered_set<T>&& s) : std::unordered_set<T>(s) {}

    Set(Set&& s) : std::unordered_set<T>((std::unordered_set<T>&&)s) {}

    Set(const std::initializer_list<T>& items) requires std::copy_constructible<T> 
        : std::unordered_set<T>(items)
    {}


    Set(const std::initializer_list<T>& items) : std::unordered_set<T>()
    {
        for (auto& item : items) add(item.clone());
    }

    Set clone() const
    {
        if constexpr (std::is_copy_constructible<T>::value)
        {
            return Set((const std::unordered_set<T>&)(*this));
        }
        else
        {
            Set<T> s;
            for (auto& item : *this) s.add(item.clone());
            return s;
        }
    }

    void operator=(Set&& s)
    {
        std::unordered_set<T>::operator=((std::unordered_set<T>&&)s);
    }

    void add(const T& val)
    {
        std::unordered_set<T>::insert(val);
    }

    bool contains(const T& val)
    {
        return std::unordered_set<T>::count(val) > 0;
    }

    void remove(const T& val)
    {
        std::unordered_set<T>::erase(val);
    }
};

template<Printable T>
std::ostream& operator<<(std::ostream& os, const Set<T>& s)
{
    int i = 0, sz = s.size();
    os << "{";
    for (const auto& it : s) {
        os << it;
        if (++i < sz) os << ", ";
    }
    os << "}";
    return os;
}

template<typename K, typename V>
class Dict final : public Object, public std::unordered_map<K, V>
{
    Dict(const Dict& d) = delete;
    void operator=(const Dict&) = delete;
    friend class List<Dict<K, V>>;
    template<typename K2, typename V2> friend class Dict;
    template<typename K2, typename V2> friend struct std::pair;
public:
    Dict() : std::unordered_map<K, V>()
    {
    }

    Dict(const std::unordered_map<K, V>& d) : std::unordered_map<K, V>(d) {}

    Dict(std::unordered_map<K, V>&& d) : std::unordered_map<K, V>(d) {}

    Dict(Dict&& d) : std::unordered_map<K, V>((std::unordered_map<K, V>&&)d) {}

    Dict(const std::initializer_list<std::pair<const K, V>>& items) 
        requires std::copy_constructible<K> && std::copy_constructible<V>
        : std::unordered_map<K, V>(items)
    {}

    Dict(const std::initializer_list<std::pair<const K, V>>& items)
        requires std::copy_constructible<K> && Clonable<V>
        : std::unordered_map<K, V>()
    {
        for (auto& item : items)
        {
            this->operator[](item.first) = item.second.clone();
        }
    }

    Dict clone() const requires std::copy_constructible<K> && std::copy_constructible<V>
    {
        return Dict((const std::unordered_map<K, V>&) (*this));
    }

    Dict clone() const requires std::copy_constructible<K> && Clonable<V>
    {
        Dict<K, V> d;
        for (const auto& it : *this)
        {
            d[it.first] = it.second.clone();
        }
        return d;
    }

    void operator=(Dict&& d)
    {
        std::unordered_map<K, V>::operator=((std::unordered_map<K, V>&&)d);
    }

    bool contains_key(const K& key)
    {
        return this->find(key) != this->end();
    }
};

template<Printable K, Printable V>
std::ostream& operator<<(std::ostream& os, const Dict<K, V>& dict) {

    os << "{";
    int i = 0, sz = dict.size();
    for (auto& it : dict) {
        os << it.first << ": " << it.second;
        if (++i < sz) os << ", ";
    }
    os << "}";

    return os;
}

class String : public Object, public std::string
{
    String(const String& s) : std::string((std::string&&)std::move(s)) {};
    void operator =(const String&) = delete;
    friend class List<String>;
    template<typename K2, typename V2> friend class Dict;
    template<typename K2, typename V2> friend struct std::pair;

public:
    String()
    {
    }

    String(const char* st) : std::string(st) {}
    String(const std::string& st) : std::string(st) {}
    String(const char ch) : std::string({ ch }) {}
    String(std::string&& st) : std::string((std::string&&)st) {}
    String(String&& st) : std::string((std::string&&)st) {}

    String clone() const
    {
        return String((std::string)(*this));
    }

    String substr(size_t start, size_t end)
    {
        if (end < start)
        {
            return String();
        }

        return String(std::string::substr(start, end - start));
    }

    template<typename T>
    void operator +=(const T& val)
    {
        if constexpr (
            std::is_same_v<T, int> || std::is_same_v<T, long> || std::is_same_v<T, long long> ||
            std::is_same_v<T, float> || std::is_same_v<T, double> ||
            std::is_same_v<T, unsigned int> || std::is_same_v<T, unsigned long> ||
            std::is_same_v<T, unsigned short> || std::is_same_v<T, short>
            )
        {
            std::string::operator+=(std::to_string(val));
        }
        else if constexpr (std::is_same_v<T, bool>)
        {
            std::string::operator+=(val ? "true" : "false");
        }
        else
        {
            std::string::operator+=(val);
        }
    }

    template<typename T>
    void append(const T& val)
    {
        operator += (val);
    }

    size_t len()
    {
        return std::string::size();
    }

    String substr(size_t start)
    {
        size_t l = len();
        if (start >= l)
        {
            return String();
        }
        return substr(start, l);
    }

    void operator = (String&& st)
    {
        std::string::operator=((std::string&&)st);
    }

    int find(const String& str)
    {
        size_t pos = std::string::find(str);
        return pos == std::string::npos ? -1 : pos;
    }

    int rfind(const String& str)
    {
        size_t pos = std::string::rfind(str);
        return pos == std::string::npos ? -1 : pos;
    }

    List<String> split(const String& delimiter)
    {
        List<String> lst;
        int pos = 0;
        String tmp = clone();
        while ((pos = tmp.find(delimiter)) >= 0) {
            lst.append(tmp.substr(0, pos));
            tmp.erase(0, pos + delimiter.length());
        }
        lst.append(std::move(tmp));
        return lst;
    }
};

template<>
struct std::hash<String>
{
    std::size_t operator()(String const& st) const noexcept
    {
        return std::hash<std::string>{}((std::string)st);
    }
};

static std::ostream& operator<<(std::ostream& os, const String& st) {
    os << "\"" << st.data() << "\"";
    return os;
}

template<typename ITR, bool ref = false>
class EnumerateIterator
{
    ITR itr;
    int index;
public:
    EnumerateIterator()
    {
        index = 0;
        itr = ITR();
    }

    EnumerateIterator(const ITR& itr) {
        this->itr = itr;
        this->index = 0;
    }

    EnumerateIterator& operator ++() {
        ++index;
        ++itr;
        return *this;
    }

    EnumerateIterator& operator += (int step)
    {
        index += step;
        itr += step;
        return *this;
    }

    bool operator != (EnumerateIterator& rhs) {
        return rhs.itr != itr;
    }

    auto operator*() const {
        if constexpr (ref)
        {
            return std::tie(index, *itr);
        }
        else
        {
            return std::make_tuple(index, *itr);
        }
    }
};


template<typename U, bool ref = false>
class Enumerate {
    Ptr<U> lst;
    using ITR = decltype(lst->begin());
public:
    Enumerate() {
        lst = (U)0;
    }

    Enumerate(U& lst) {
        this->lst = &lst;
    }

    EnumerateIterator<ITR, ref> begin() {
        return EnumerateIterator<ITR, ref>(lst->begin());
    }

    EnumerateIterator<ITR, ref> end() {
        return EnumerateIterator<ITR, ref>(lst->end());
    }

    size_t size() const
    {
        return lst->size();
    }

    auto par_iter(int n_thread)
    {
        return get_par_iter(*this, n_thread);
    }
};

template<typename U>
auto enumerate(U& lst) {
    return Enumerate(lst);
}

template<typename U>
auto enumerate_ref(U& lst) {
    return Enumerate<U, true>(lst);
}

template<typename ITR1, typename ITR2, bool ref = false>
class Zip2Iterator
{
    ITR1 itr1;
    ITR2 itr2;
public:
    Zip2Iterator()
    {
    }

    Zip2Iterator(const ITR1& itr1, const ITR2& itr2)
    {
        this->itr1 = itr1;
        this->itr2 = itr2;
    }

    Zip2Iterator& operator ++()
    {
        ++itr1;
        ++itr2;
        return *this;
    }

    Zip2Iterator& operator += (int step)
    {
        itr1 += step;
        itr2 += step;
        return *this;
    }

    bool operator != (Zip2Iterator& rhs) {
        return rhs.itr1 != itr1 && rhs.itr2 != itr2;
    }

    auto operator*() const {
        if constexpr (ref)
        {
            return std::forward_as_tuple(*itr1, *itr2);
        }
        else
        {
            return std::make_tuple(*itr1, *itr2);
        }
    }
};

template<typename U1, typename U2, bool ref = false>
class Zip2 {
    Ptr<U1> lst1;
    Ptr<U2> lst2;
    using ITR1 = decltype(lst1->begin());
    using ITR2 = decltype(lst2->begin());
public:
    Zip2()
    {
    }

    Zip2(U1& lst1, U2& lst2) {
        this->lst1 = &lst1;
        this->lst2 = &lst2;
    }

    Zip2Iterator<ITR1, ITR2, ref> begin() {
        return Zip2Iterator<ITR1, ITR2, ref>(lst1->begin(), lst2->begin());
    }

    Zip2Iterator<ITR1, ITR2, ref> end() {
        return Zip2Iterator<ITR1, ITR2, ref>(lst1->end(), lst2->end());
    }


    size_t size() const
    {
        size_t sz1 = lst1->size();
        size_t sz2 = lst2->size();
        return (sz1 < sz2) ? sz1 : sz2;
    }

    auto par_iter(int n_thread)
    {
        return get_par_iter(*this, n_thread);
    }
};

template<typename U1, typename U2>
auto zip(U1& lst1, U2& lst2) {
    return Zip2(lst1, lst2);
}

template<typename U1, typename U2>
auto zip_ref(U1& lst1, U2& lst2) {
    return Zip2<U1, U2, true>(lst1, lst2);
}

template<typename ITR1, typename ITRM, bool ref = false>
class ZipMultiIterator
{
    ITR1 itr1;
    ITRM itrm;
public:
    ZipMultiIterator()
    {
    }

    ZipMultiIterator(const ITR1& itr1, const ITRM& itrm)
    {
        this->itr1 = itr1;
        this->itrm = itrm;
    }

    ZipMultiIterator& operator ++()
    {
        ++itr1;
        ++itrm;
        return *this;
    }

    ZipMultiIterator& operator += (int step)
    {
        itr1 += step;
        itrm += step;
        return *this;
    }

    bool operator != (ZipMultiIterator& rhs)
    {
        return itr1 != rhs.itr1 && itrm != rhs.itrm;
    }

    auto operator*() const {
        if constexpr (ref)
        {
            return std::tuple_cat(std::tie(*itr1), *itrm);
        }
        else
        {
            return std::tuple_cat(std::make_tuple(*itr1), *itrm);
        }
    }
};

template<bool ref>
class Phantom {};

template<typename U, typename Z, bool ref = false>
class ZipMulti {
    Ptr<U> lst;
    Z z;
    using ITR1 = decltype(lst->begin());
    using ITRM = decltype(z.begin());
public:
    ZipMulti() {}

    ZipMulti(U& lst, const Z& z, const Phantom<ref>& phantom) {
        this->lst = &lst;
        this->z = z;
    }

    ZipMultiIterator<ITR1, ITRM, ref> begin() {
        return ZipMultiIterator<ITR1, ITRM, ref>(lst->begin(), z.begin());
    }

    ZipMultiIterator<ITR1, ITRM, ref> end() {
        return ZipMultiIterator<ITR1, ITRM, ref>(lst->end(), z.end());
    }

    size_t size() const
    {
        size_t sz1 = lst->size();
        size_t sz2 = z.size();
        return (sz1 < sz2) ? sz1 : sz2;
    }

    auto par_iter(int n_thread)
    {
        return get_par_iter(*this, n_thread);
    }
};

template<typename U, typename... Args>
auto zip(U& lst, Args&... args)
{
    return ZipMulti(lst, zip(args...), Phantom<false>());
}

template<typename U, typename... Args>
auto zip_ref(U& lst, Args&... args)
{
    return ZipMulti(lst, zip_ref(args...), Phantom<true>());
}

#endif
