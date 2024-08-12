#ifndef __ARRAY_H__
#define __ARRAY_H__

#include "base.h"
#include <math.h>
#include <typeinfo>

constexpr size_t NUM_OP = 15;
constexpr size_t OP_ADD = 1;
constexpr size_t OP_SUB = 2;
constexpr size_t OP_MUL = 3;
constexpr size_t OP_DIV = 4;
constexpr size_t OP_MOD = 5;
constexpr size_t OP_POW = 6;
constexpr size_t OP_AND = 7;
constexpr size_t OP_OR  = 8;
constexpr size_t OP_XOR = 9;
constexpr size_t OP_GT  = 10;
constexpr size_t OP_GTE = 11;
constexpr size_t OP_LT  = 12;
constexpr size_t OP_LTE = 13;
constexpr size_t OP_EQ  = 14;
constexpr size_t OP_NE  = 15;

constexpr size_t NUM_REDUCE = 12;
constexpr size_t REDUCE_MAX = 1;
constexpr size_t REDUCE_MIN = 2;
constexpr size_t REDUCE_SUM = 3;
constexpr size_t REDUCE_AND = 4;
constexpr size_t REDUCE_OR  = 5;
constexpr size_t REDUCE_XOR = 6;
constexpr size_t REDUCE_MEAN = 7;
constexpr size_t REDUCE_MEDIAN = 8;
constexpr size_t REDUCE_STD = 9;
constexpr size_t REDUCE_ARGMIN = 10;
constexpr size_t REDUCE_ARGMAX = 11;
constexpr size_t REDUCE_PERCENTILE = 12;

constexpr size_t NUM_FUNC = 16;
constexpr size_t FUNC_ABS = 1;
constexpr size_t FUNC_ROUND = 2;
constexpr size_t FUNC_CEIL = 3;
constexpr size_t FUNC_FLOOR = 4;
constexpr size_t FUNC_SQUARE = 5;
constexpr size_t FUNC_SQRT = 6;
constexpr size_t FUNC_EXP = 7;
constexpr size_t FUNC_LOG = 8;
constexpr size_t FUNC_SIN = 9;
constexpr size_t FUNC_COS = 10;
constexpr size_t FUNC_TAN = 11;
constexpr size_t FUNC_ASIN = 12;
constexpr size_t FUNC_ACOS = 13;
constexpr size_t FUNC_ATAN = 14;
constexpr size_t FUNC_NEG = 15;
constexpr size_t FUNC_NOT = 16;

template<typename T>
static void* _typeId() //this function is instantiated for every different type
{
    static T* marker = NULL; //thus this static variable will be created for each TypeIdNoRTTI<T> separately
    return &marker; //it's address is unique identifier of TypeIdNoRTTI<T> type
}   

namespace np {
    const int UINT8 = 1;
    const int INT8 = 2;
    const int UINT16 = 3;
    const int INT16 = 4;
    const int UINT32 = 5;
    const int INT32 = 6;
    const int UINT64 = 7;
    const int INT64 = 8;
    const int BOOL = 9;
    const int FLOAT32 = 10;
    const int FLOAT64 = 11;
};

template<typename T>
int getTypeId()
{
    auto& type_id = typeid(T);
    if(type_id == typeid(uint8_t)) return np::UINT8;
    if(type_id == typeid(int8_t)) return np::INT8;
    if(type_id == typeid(uint16_t)) return np::UINT16;
    if(type_id == typeid(int16_t)) return np::INT16;
    if(type_id == typeid(uint32_t)) return np::UINT32;
    if(type_id == typeid(int32_t)) return np::INT32;
    if(type_id == typeid(uint64_t)) return np::UINT64;
    if(type_id == typeid(int64_t)) return np::INT64;
    if(type_id == typeid(bool)) return np::BOOL;
    if(type_id == typeid(float)) return np::FLOAT32;
    if(type_id == typeid(double)) return np::FLOAT64;
}

template<typename T>
void swap(T* a, T* b)
{
    T t = *a;
    *a = *b;
    *b = t;
}

template<typename T>
void swap_with_idx(T* a, T* b, int* ia, int* ib)
{
    T t = *a;
    *a = *b;
    *b = t;
    int i = *ia;
    *ia = *ib;
    *ib = i;
}

template<typename T>
int partition(T arr[], int step, int indexes[], int low, int high)
{
    // Choosing the pivot
    T pivot = arr[high*step];
 
    // Index of smaller element and indicates
    // the right position of pivot found so far
    int i = (low - 1);
 
    for (int j = low; j <= high - 1; j++) {
 
        // If current element is smaller than the pivot
        if (arr[j*step] < pivot) {
 
            // Increment index of smaller element
            i++;
            if(indexes)
            {
                swap_with_idx(&arr[i*step], &arr[j*step], &indexes[i], &indexes[j]);
            }else
            {
                swap(&arr[i*step], &arr[j*step]);
            }
        }
    }
    if(indexes)
    {
        swap_with_idx(&arr[(i+1)*step], &arr[high*step], &indexes[i+1], &indexes[high]);
    }else
    {
        swap(&arr[(i+1)*step], &arr[high*step]);
    }
    return (i + 1);
}

template<typename T>
void quickSort(T arr[], int step, int indexes[], int low, int high)
{
    if (low < high) {
 
        // pi is partitioning index, arr[p]
        // is now at right place
        int pi = partition(arr, step, indexes, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(arr, step, indexes, low, pi - 1);
        quickSort(arr, step, indexes, pi + 1, high);
    }
}

struct Slice 
{
    Integer start = Integer::none();
    Integer stop = Integer::none();
    Integer step = Integer::none();

    Slice()
    {

    }

    Slice(Integer start, Integer stop, Integer step=Integer::none())
    {
        this->start = start;
        this->stop = stop;
        this->step = step;
    }

    static Slice all()
    {
        return Slice();
    }

    static Slice head(int stop)
    {
        return Slice(Integer::none(), Integer(stop), Integer::none());
    }

    static Slice head(int stop, int step)
    {
        return Slice(Integer::none(), Integer(stop), Integer(step));
    }

    static Slice tail(int start)
    {
        return Slice(Integer(start), Integer::none(), Integer::none());
    }

    static Slice tail(int start, int step)
    {
        return Slice(Integer(start), Integer::none(), Integer(step));
    }
};

static Slice resolve_index(Slice idx, int size) 
{
    int step = idx.step.is_null ? 1 : idx.step.value;
    int start = idx.start.is_null? (step > 0? 0 : size-1) : idx.start.value;
    int stop = idx.stop.is_null? (step > 0? size : -1) : idx.stop.value;

    if(start < 0) start += size;
    if(stop < 0 && !idx.stop.is_null) stop += size;
    if(start < 0) start = 0;
    if(stop > size) stop = size;
    
    if((step > 0 && stop < start) || (step < 0 && stop > start)){
        stop = start;
    }

    if(step > 0){
        stop = start + ((stop - start - 1) / step) * step + step;
    }else if(step < 0){
        stop = start - ((start - stop - 1) / (-step)) * (-step) + step;
    }
    return Slice(start, stop, step);
}

template <typename... T>
struct Tuple;

template<typename T1, typename T2>
struct Tuple<T1, T2>
{
    T1 item1;
    T2 item2;

    Tuple(T1 item1, T2 item2)
    {
        this->item1 = item1;
        this->item2 = item2;
    }

    constexpr int dim() const 
    {
        return 2;
    }

    int size() const
    {
        return item1 * item2;
    }
};

template<typename T1, typename T2, typename  T3>
struct Tuple<T1, T2, T3>
{
    T1 item1;
    T2 item2;
    T3 item3;

    Tuple(T1 item1, T2 item2, T3 item3)
    {
        this->item1 = item1;
        this->item2 = item2;
        this->item3 = item3;
    }

    constexpr int dim() const 
    {
        return 3;
    }

    int size() const
    {
        return item1 * item2 * item3;
    }
};

template<typename T1, typename T2, typename  T3, typename  T4>
struct Tuple<T1, T2, T3, T4>
{
    T1 item1;
    T2 item2;
    T3 item3;
    T4 item4;

    Tuple(T1 item1, T2 item2, T3 item3, T4 item4)
    {
        this->item1 = item1;
        this->item2 = item2;
        this->item3 = item3;
        this->item4 = item4;
    }

    constexpr int dim() const 
    {
        return 4;
    }

    int size() const
    {
        return item1 * item2 * item3 * item4;
    }
};

template<typename T>
class Array;

template<typename T>
class Array2;

template<typename T>
class Array3;

template<typename T>
class Array4;

template<typename T>
class Array2Iterator;


template<typename T, bool compact=0>
class ArrayIterator {
    int _step = 0;
    T* _data = NULL;
public:
    ArrayIterator()
    {
    }
    ArrayIterator(T* data, int step){
        _step = step;
        _data = data;
        if constexpr(compact) step = 1;
    }

    ArrayIterator& operator++(){
        if constexpr(compact) _data ++;
        else _data += _step;
        return *this;
    }
    
    ArrayIterator& operator+=(int step){
        _data += _step * step;
        return *this;
    }

    bool operator != (ArrayIterator& rhs) {
        return _data != rhs._data;
    }

    T& operator*() const {
        return *_data;
    } 
};

template<typename T, int func>
inline static T _unary_ops(const T& val)
{
    if constexpr(func == FUNC_ABS) return ::fabs(val);
    if constexpr(func == FUNC_ROUND) return ::round(val);
    if constexpr(func == FUNC_CEIL) return ::ceil(val);
    if constexpr(func == FUNC_FLOOR) return ::floor(val);
    if constexpr(func == FUNC_SQUARE) return val*val;
    if constexpr(func == FUNC_SQRT) return ::sqrt(val);
    if constexpr(func == FUNC_EXP) return ::exp(val);
    if constexpr(func == FUNC_LOG) return ::log(val);
    if constexpr(func == FUNC_SIN) return ::sin(val);
    if constexpr(func == FUNC_COS) return ::cos(val);
    if constexpr(func == FUNC_TAN) return ::tan(val);
    if constexpr(func == FUNC_ASIN) return ::asin(val);
    if constexpr(func == FUNC_ACOS) return ::acos(val);
    if constexpr(func == FUNC_ATAN) return ::atan(val);
    if constexpr(func == FUNC_NEG) return -val;
    if constexpr(func == FUNC_NOT) return ~val;
}

template<typename T, typename U, int ops>
inline static U _bin_ops(const T& l, const T& r)
{
    if constexpr(ops == OP_ADD) return l + r;
    if constexpr(ops == OP_SUB) return l - r;
    if constexpr(ops == OP_MUL) return l * r;
    if constexpr(ops == OP_DIV) return l / r;
    if constexpr(ops == OP_MOD) return l % r;
    if constexpr(ops == OP_POW) return ::pow(l, r);
    if constexpr(ops == OP_AND) return l & r;
    if constexpr(ops == OP_OR)  return l | r;
    if constexpr(ops == OP_XOR) return l ^ r;
    if constexpr(ops == OP_GT)  return l > r;
    if constexpr(ops == OP_GTE) return l >= r;
    if constexpr(ops == OP_LT)  return l < r;
    if constexpr(ops == OP_LTE) return l <= r;
    if constexpr(ops == OP_EQ)  return l == r;
    if constexpr(ops == OP_NE)  return l != r;
}

template<typename T, int ops>
inline static void _bin_ops_assign(T& l, const T& r) 
{
    if constexpr(ops == OP_ADD) l += r;
    if constexpr(ops == OP_SUB) l -= r;
    if constexpr(ops == OP_MUL) l *= r;
    if constexpr(ops == OP_DIV) l /= r;
    if constexpr(ops == OP_MOD) l %= r;
    if constexpr(ops == OP_AND) l &= r;
    if constexpr(ops == OP_OR)  l |= r;
    if constexpr(ops == OP_XOR) l ^= r;

}

struct CArray 
{
    void* data;
    isize size;
    isize step;
    isize dtype;
};

template<typename BOOL, typename INT>
Array<INT> mask_to_indexes(const Array<BOOL>& masks);


template<typename T>
class Array : public Object 
{
    template<typename U> friend class Array;
    template<typename U> friend class Array2;
    template<typename U> friend class Array3;
    template<typename U> friend class Array4;

    friend class Array2Iterator<T>;

    int _size = 0;
    int _step = 0;
    T* _data = NULL;
    bool _is_view = false;
    const size_t* _p_signature = NULL;

    void operator = (const Array&) = delete;

    Array(const Array&) = delete;

    const size_t* get_p_signature() const
    {
        return _is_view? _p_signature: &_signature;
    }

    void _check_valid() const 
    {
        if(_is_view)
        {
             #ifdef DEBUG
            if(*get_p_signature() != get_signature() || !_data){
                panic("Array is invalid\n");
            }
            #endif
        }
    }

    Array(T* data, int size, int step, bool is_view=false, const size_t* p_signature=NULL)
    {
        _data = data;
        _size = size;
        _step = step;
        _is_view = is_view;
        _p_signature = p_signature;
        if(p_signature)
        {
            _signature = *p_signature;
        }
    }
    
public:
    Array()
    {
    }

    Array(const CArray* carr, char* p_signature): Array((T*) carr->data, carr->size, carr->step, true, p_signature)
    {

    }

    Array(Array&& oth)
    {
        _data = oth._data;
        _step = oth._step;
        _size = oth._size;
        _signature = oth._signature;
        _p_signature = oth._p_signature;
        _is_view = oth._is_view;
        oth._data = NULL;
        oth._size = 0;
    }

    void operator = (Array&& oth)
    {
        _check_valid();
        oth._check_valid();

        if(!_is_view && _data) delete[] _data;
        _data = oth._data;
        _step = oth._step;
        _size = oth._size;
        _signature = oth._signature;
        _p_signature = oth._p_signature;
        _is_view = oth._is_view;
        oth._data = NULL;
        oth._size = 0;
    }

    static Array new_array(T* raw, int n)
    {
        T* data = new T[n];
        for(int i = 0; i < n; i++) data[i] = raw[i];
        return Array(data, n, 1);
    }

    static Array from_raw(T* raw, int n)
    {
        return Array(raw, n, 1);
    }

    static Array full(int n, T value)
    {
        T* data = new T[n];
        for(int i = 0; i < n; i++) data[i] = value;
        return Array(data, n, 1);
    }

    static Array zeros(int n)
    {
        return full(n, (T) 0);
    }

    static Array ones(int n)
    {
        return full(n, (T) 1);
    }

    static Array from_list(const List<T>& lst)
    {
        int n = lst.size();
        T* data = new T[n];
        for(int i = 0; i < n; i++) data[i] = lst[i];
        return Array(data, n, 1);
    }

    template<typename INT>
    void get_by_indexes_to(const Array<INT>& indexes, Array<T>& out)
    {
        _check_valid();

        for(int i = 0; i < out.size(); i++) 
        {
            out._at(i) = operator[](indexes._at(i));
        }
    }

    Array<T> operator[](const Array<int>& indexes)
    {
        auto res = Array<T>::zeros(indexes.size());
        get_by_indexes_to(indexes, res);
        return move(res);
    }

    template<typename INT>
    void set_by_indexes(const Array<INT>& indexes, const T& val)
    {
        _check_valid();

        for(INT idx: indexes) 
        {
            operator[](idx) = val;
        }
    }

    template<typename INT>
    void set_by_indexes(const Array<INT>& indexes, const Array<T>& val)
    {
        _check_valid();

        assert(indexes.size() == val.size());

        for(int i = 0; i < indexes.size(); i++) 
        {
            operator[](indexes._at(i)) = val._at(i);
        }
    }

    template<typename DATA>
    void set(const Array<int>& indexes, const DATA& val)
    {
        set_by_indexes(indexes, val);
    }

    Array<T> operator[](const Array<bool>& masks)
    {
        _check_valid();
        assert(masks.size() == _size);
        return operator[](mask_to_indexes<bool, int>(masks));
    }

    template<typename DATA>
    void set(const Array<bool>& masks, const DATA& val) 
    {
        _check_valid();
        assert(masks.size() == _size);
        set(mask_to_indexes<bool, int>(masks), val);
    }

    Array2<T> reshape(Tuple<int, int> new_shape) 
    {
        _check_valid();
        if(new_shape.item1 * new_shape.item2 != _size) 
        {
            panic("Cannot reshape array of size %d to new shape (%d,%d)\n", 
                    _size, new_shape.item1, new_shape.item2);
        }
        Array2<T> res = Array2<T>::zeros(new_shape);
        if(is_compact())
        {
            for(int i = 0; i < _size; i++)
            {
                res._data[i] = _at_compact(i);
            }
        }
        else
        {
            for(int i = 0; i < _size; i++)
            {
                res._data[i] = _at(i);
            }
        }
        return res;
    }

    const T* get_ptr() const
    {
        _check_valid();
        return _data;
    }

    T* get_ptr()
    {
        _check_valid();
        return _data;
    }

    int size() const 
    {
        _check_valid();
        return _size;
    }

    int step() const 
    {
        _check_valid();
        return _step;
    }

    bool is_view() const 
    {
        return _is_view;
    }

    bool is_compact() const
    {
        return _step == 1;
    }

    T& _at(int idx) const
    {
        return _data[_step * idx];
    }

    T& _at_compact(int idx) const
    {
        return _data[idx];
    }

    T& operator[](int idx_) const
    {
        _check_valid();

        int idx = idx_ < 0 ? _size + idx_ : idx_;
        #ifdef DEBUG
        if(idx < 0 || idx >= _size)
        {
            panic("Index %d is out of range %d\n", idx_, _size);
        }
        #endif
        return _at(idx);
    }

    template<typename U>
    Array<U> as_type() const {
        _check_valid();

        U* data = new U[_size];
        if(is_compact())
        {
            for(int i = 0; i < _size; i++) data[i] = (U) _at_compact(i);
        }
        else
        {
            for(int i = 0; i < _size; i++) data[i] = (U) _at(i);
        }
        return Array<U>(data, _size, 1);
    }

    Array<T> copy() const
    {
        return as_type<T>();
    }

    T* into_raw() 
    {
        _check_valid();
        T* data = _data;
        _data = NULL;
        _size = 0;
        return data;
    }

    template<typename U>
    void copy_from(const Array<U>& oth)
    {
        _check_valid();
        oth._check_valid();

        if(_size != oth._size)
        {
            panic("Array size mismatch %d vs %d \n", _size, oth._size);
        }

        if(is_compact() && oth.is_compact())
        {
            for(int i = 0; i < _size; i++)
            {
                _at_compact(i) = (T) oth._at_compact(i);
            }
        }
        else
        {
            for(int i = 0; i < _size; i++)
            {
                _at(i) = (T) oth._at(i);
            }
        }
    }

    template<typename U>
    void assign(const Array<U>& oth)
    {
        copy_from(oth);
    }

    void set_all(const T& value)
    {
        _check_valid();
        if(is_compact())
        {
            for(int i = 0; i < _size; i++) _at_compact(i) = value;
        }
        else
        {
            for(int i = 0; i < _size; i++) _at(i) = value;
        }
    }

    void assign(const T& value)
    {
        set_all(value);
    }

    ArrayIterator<T> begin() const
    {
        _check_valid();
        return ArrayIterator<T>(_data, _step);
    }

    ArrayIterator<T> end() const
    {
        _check_valid();
        return ArrayIterator<T>(_data + _size * _step, _step);
    }

    Array<T> view() const
    {
        return Array<T>(_data, _size, _step, true, get_p_signature());
    }

    Array<T> slice(Integer start=0, Integer stop=Integer::none(), Integer step=Integer::none())
    {
        return get(Slice(start, stop, step));
    }

    Array<T> get(Slice idx)
    {
        return operator[](idx);
    }

    Array<T> operator[](Slice idx)
    {
        _check_valid();

        idx = resolve_index(idx, _size);
        int start = idx.start.value;
        int stop =  idx.stop.value;
        int step =  idx.step.value;
    
        assert(step != 0);

        return Array<T>(
            _data + _step * start,          // _data
            (stop - start)/step,            // _size
            _step * step,                   // _step
            true,
            get_p_signature()
        );

    }

    void set(Slice idx, const Array<T>& val)
    {
        operator[](idx).copy_from(val);
    }

    void set(Slice idx, const T& val)
    {
        operator[](idx).set_all(val);
    }

    template<int func>
    void unary_ops_compute(Array<T>& res) const
    {
        _check_valid();
        res._check_valid();

        assert(func <= NUM_FUNC);

        if(is_compact() && res.is_compact())
        {
            for(int i = 0; i < _size; i++)
            {
                res._at_compact(i) = _unary_ops<T, func>(_at_compact(i));
            }
        }
        else
        {
            for(int i = 0; i < _size; i++)
            {
                res._at(i) = _unary_ops<T, func>(_at(i));
            }
        }
    }

    template<int func>
    Array<T> unary_ops() const
    {
        auto res = Array<T>::zeros(_size);
        unary_ops_compute<func>(res);
        return res;
    }

    Array<T> operator ~() 
    {
        return unary_ops<FUNC_NOT>();
    }

    Array<T> operator -() 
    {
        return unary_ops<FUNC_NEG>();
    }

    Array<T> abs() const {
        return unary_ops<FUNC_ABS>();
    }

    Array<T> round() const {
        return unary_ops<FUNC_ROUND>();
    }

    Array<T> ceil() const {
        return unary_ops<FUNC_CEIL>();
    }

    Array<T> floor() const {
        return unary_ops<FUNC_FLOOR>();
    }

    Array<T> square() const {
        return unary_ops<FUNC_SQUARE>();
    }

    Array<T> sqrt() const {
        return unary_ops<FUNC_SQRT>();
    }

    Array<T> exp() const {
        return unary_ops<FUNC_EXP>();
    }

    Array<T> log() const {
        return unary_ops<FUNC_LOG>();
    }

    Array<T> sin() const {
        return unary_ops<FUNC_SIN>();
    }

    Array<T> cos() const {
        return unary_ops<FUNC_COS>();
    }

    Array<T> tan() const {
        return unary_ops<FUNC_TAN>();
    }

    Array<T> asin() const {
        return unary_ops<FUNC_ASIN>();
    }

    Array<T> acos() const {
        return unary_ops<FUNC_ACOS>();
    }

    Array<T> atan() const {
        return unary_ops<FUNC_ATAN>();
    }

    template<int ops>
    void bin_ops_assign(const T& rhs) {
        _check_valid();

        if constexpr(ops > NUM_OP)
        {
            panic("Invalid operator\n");
        }

        if(is_compact())
        {
            for(int i = 0; i < _size; i++)
            {
                _bin_ops_assign<T, ops>(_at_compact(i), rhs);
            }
        }
        else
        {
            int i = 0;
            for(; i+3 < _size; i+=4)
            {
                _bin_ops_assign<T, ops>(_at(i), rhs);
                _bin_ops_assign<T, ops>(_at(i+1), rhs);
                _bin_ops_assign<T, ops>(_at(i+2), rhs);
                _bin_ops_assign<T, ops>(_at(i+3), rhs);
            }
            for(; i < _size; i++)
            {
                _bin_ops_assign<T, ops>(_at(i), rhs);
            }
        }
        
    }

    void operator += (const T&rhs) {
        bin_ops_assign<OP_ADD>(rhs);
    }

    void operator -=(const T& rhs) {
        bin_ops_assign<OP_SUB>(rhs);
    }

    void operator *=(const T& rhs) {
        bin_ops_assign<OP_MUL>(rhs);
    }

    void operator /=(const T& rhs) {
        bin_ops_assign<OP_DIV>(rhs);
    }

    void operator %=(const T& rhs) {
        bin_ops_assign<OP_MOD>(rhs);
    }

    template<int ops>
    void bin_ops_assign(const Array<T>& rhs) {
        _check_valid();
        rhs._check_valid();

        if constexpr(ops > NUM_OP)
        {
            panic("Invalid operator : %d\n", ops);
        }

        if(rhs._size != _size && rhs._size != 1 )
        {
            panic("Dimensions mismatch: %d vs %d\n", _size, rhs._size);
        }

        if(rhs._size == _size)
        {
            if(is_compact() && rhs.is_compact())
            {
                for(int i = 0; i < _size; i++)
                {
                    _bin_ops_assign<T, ops>(_at_compact(i), rhs._at_compact(i));
                }
            }
            else
            {
                int i = 0;
                for(; i+3 < _size; i+=4)
                {
                    _bin_ops_assign<T, ops>(_at(i), rhs._at(i));
                    _bin_ops_assign<T, ops>(_at(i+1), rhs._at(i+1));
                    _bin_ops_assign<T, ops>(_at(i+2), rhs._at(i+2));
                    _bin_ops_assign<T, ops>(_at(i+3), rhs._at(i+3));
                }
                for(;i < _size; i++)
                {
                    _bin_ops_assign<T, ops>(_at(i), rhs._at(i));
                }
            }
        }
        else
        {
            bin_ops_assign<ops>(rhs._at(0));
        }
    }

    void operator +=(const Array<T>& rhs) {
        bin_ops_assign<OP_ADD>(rhs);
    }

    void operator -=(const Array<T>& rhs) {
        bin_ops_assign<OP_SUB>(rhs);
    }

    void operator *=(const Array<T>& rhs) {
        bin_ops_assign<OP_MUL>(rhs);
    }

    void operator /=(const Array<T>& rhs) {
        bin_ops_assign<OP_DIV>(rhs);
    }

    void operator %=(const Array<T>& rhs) {
        bin_ops_assign<OP_MOD>(rhs);
    }

    void operator &=(const Array<T>& rhs) {
        bin_ops_assign<OP_AND>(rhs);
    }

    void operator |=(const Array<T>& rhs) {
        bin_ops_assign<OP_OR>(rhs);
    }

    void operator ^=(const Array<T>& rhs) {
        bin_ops_assign<OP_XOR>(rhs);
    }

    template<typename U, int ops>
    static void bin_ops_compute(const Array<T>& lhs, const T& rhs, Array<U>& result) 
    {
        if(result.is_compact() && lhs.is_compact())
        {
            for(int i = 0;i < result._size; i++)
            {
                result._at_compact(i) = _bin_ops<T, U, ops>(lhs._at_compact(i), rhs);
            }
        }
        else
        {
            for(int i = 0; i < result._size; i++)
            {
                result._at(i) = _bin_ops<T, U, ops>(lhs._at(i), rhs);
            }
        }
    }

    template<typename U, int ops>
    Array<U> bin_ops(const T& rhs) const {
        _check_valid();

        if constexpr(ops > NUM_OP)
        {
            panic("Invalid operator: %d\n", ops);
        }

        auto arr = Array<U>::zeros(_size);

        bin_ops_compute<U, ops>(*this, rhs, arr);

        return arr;
    }

    Array<T> operator +(const T& rhs) {
        return bin_ops<T,OP_ADD>(rhs);
    }

    Array<T> operator -(const T& rhs) {
        return bin_ops<T,OP_SUB>(rhs);
    }

    Array<T> operator *(const T& rhs) {
        return bin_ops<T,OP_MUL>(rhs);
    }

    Array<T> operator /(const T& rhs) {
        return bin_ops<T,OP_DIV>(rhs);
    }

    Array<T> operator %(const T& rhs) {
        return bin_ops<T,OP_MOD>(rhs);
    }

    Array<T> pow(const T& rhs) const {
        return bin_ops<T, OP_POW>(rhs); 
    }

    Array<T> operator &(const T& rhs) {
        return bin_ops<T,OP_AND>(rhs);
    }

    Array<T> operator |(const T& rhs) {
        return bin_ops<T,OP_OR>(rhs);
    }

    Array<T> operator ^(const T& rhs) {
        return bin_ops<T,OP_XOR>(rhs);
    }

    Array<bool> operator >(const T& rhs) {
        return bin_ops<bool,OP_GT>(rhs);
    }

    Array<bool> operator >=(const T& rhs) {
        return bin_ops<bool,OP_GTE>(rhs);
    }

    Array<bool> operator <(const T& rhs) {
        return bin_ops<bool,OP_LT>(rhs);
    }

    Array<bool> operator <=(const T& rhs) {
        return bin_ops<bool,OP_LTE>(rhs);
    }

    Array<bool> operator ==(const T& rhs) {
        return bin_ops<bool,OP_EQ>(rhs);
    }

    Array<bool> operator !=(const T& rhs) {
        return bin_ops<bool,OP_NE>(rhs);
    }

    template<typename U, int ops>
    static void bin_ops_compute(const Array<T>& lhs, const Array<T>& rhs, 
            Array<U>& result) {

        if(lhs._size == 1)
        {
            T l = lhs._at(0);
            if(result.is_compact() && rhs.is_compact())
            {
                for(int i = 0; i < result._size; i++)
                {
                    result._at_compact(i) = _bin_ops<T, U, ops>(l, rhs._at_compact(i));
                }
            }
            else
            {
                for(int i = 0; i < result._size;i++)
                {
                    result._at(i) = _bin_ops<T, U, ops>(l, rhs._at(i));
                }
            }
        }
        else if(rhs._size == 1)
        {
            bin_ops_compute<U, ops>(lhs, rhs._at(9), result);
        }
        else
        {
            if(result.is_compact() && lhs.is_compact() && rhs.is_compact())
            {
                for(int i = 0; i < result._size; i++)
                {
                    result._at_compact(i) = _bin_ops<T, U, ops>(lhs._at_compact(i), rhs._at_compact(i));
                }
            }
            else
            {
                for(int i =0; i < result._size; i++)
                {
                    result._at(i) = _bin_ops<T, U, ops>(lhs._at(i), rhs._at(i));
                }
            }
        }
    }

    template<typename U, int ops>
    Array<U> bin_ops(const Array<T>& rhs) const
    {
        if constexpr(ops > NUM_OP)
        {
            panic("Invalid operator\n");
        }

        if(_size != 1 && rhs._size != 1 && _size != rhs._size)
        {
            panic("Dimensions mismatch\n");
        }

        int sz = _size > rhs._size? _size : rhs._size;
        Array<U> arr = Array<U>::zeros(sz);
        bin_ops_compute<U, ops>(*this, rhs, arr);
        return arr;
    }

    Array<T> operator +(const Array<T>& rhs) const {
        return bin_ops<T, OP_ADD>(rhs);
    }

    Array<T> operator -(const Array<T>& rhs) const {
        return bin_ops<T, OP_SUB>(rhs);
    }

    Array<T> operator *(const Array<T>& rhs) const {
        return bin_ops<T, OP_MUL>(rhs);
    }

    Array<T> operator /(const Array<T>& rhs) const {
        return bin_ops<T, OP_DIV>(rhs);
    }

    Array<T> operator %(const Array<T>& rhs) const {
        return bin_ops<T, OP_MOD>(rhs);
    }

    Array<T> pow(const Array<T>& rhs) const {
        return bin_ops<T, OP_POW>(rhs); 
    }

    Array<T> operator &(const Array<T>& rhs)  const{
        return bin_ops<T, OP_AND>(rhs);
    }

    Array<T> operator |(const Array<T>& rhs) const{
        return bin_ops<T, OP_OR>(rhs);
    }

    Array<T> operator ^(const Array<T>& rhs) const{
        return bin_ops<T, OP_XOR>(rhs);
    }

    Array<bool> operator >(const Array<T>& rhs) const{
        return bin_ops<bool, OP_GT>(rhs);
    }

    Array<bool> operator >=(const Array<T>& rhs) const{
        return bin_ops<bool, OP_GTE>(rhs);
    }

    Array<bool> operator <(const Array<T>& rhs) const{
        return bin_ops<bool, OP_LT>(rhs);
    }

    Array<bool> operator <=(const Array<T>& rhs) const{
        return bin_ops<bool, OP_LTE>(rhs);
    }

    Array<bool> operator ==(const Array<T>& rhs) const{
        return bin_ops<bool, OP_EQ>(rhs);
    }

    Array<bool> operator !=(const Array<T>& rhs) const{
        return bin_ops<bool, OP_NE>(rhs);
    }

    double percentile(double p) const
    {
        if(_size == 0)
        {
            panic("Percentile of empty array\n");
        }
        auto arr = copy();
        arr.sort();
        double index = (p/100) * (_size-1);
        int index_int = (int) index;
        double delta = index - index_int;
        if(index_int == _size-1)
        {
            return _at(index_int);
        }
        return (1.0 - delta) * arr._at(index_int) + delta * arr._at(index_int+1);
    }

    double reduce_median() const 
    {
        return percentile(50.0);
    }

    double reduce_mean() const 
    {
        if(_size == 0)
        {
            panic("Mean of empty array\n");
        }
        return reduce_ops<double, REDUCE_SUM>() / _size;
    }

    template<typename U, int ops>
    U reduce_ops() const {
        if constexpr(ops > NUM_REDUCE)
        {
            panic("Invalid reduce operator\n");
        }

        if constexpr(ops == REDUCE_MAX || ops == REDUCE_MIN || ops == REDUCE_AND || ops == REDUCE_OR || ops == REDUCE_XOR) {
            T val = _at(0);
            for(int i = 0; i < _size; i++){
                T tmp = _at(i);
                if constexpr(ops == REDUCE_MAX) {
                    if(tmp > val) val = tmp;
                }
                if constexpr(ops == REDUCE_MIN) {
                    if(tmp < val) val = tmp;
                }
                if constexpr(ops == REDUCE_AND) {
                    tmp &= val;
                }
                if constexpr(ops == REDUCE_OR) {
                    tmp |= val;
                }
                if constexpr(ops == REDUCE_XOR) {
                    tmp ^= val;
                }
            }
            return val;
        }

        if constexpr(ops == REDUCE_SUM ) {
            U val = U(0);
            for(int i = 0; i < _size; i++){
                val += _at(i);
            }
            return val;
        }

        if constexpr(ops == REDUCE_MEAN) 
        {
            return (U) reduce_mean();
        }

        if constexpr(ops == REDUCE_MEDIAN) 
        {
            return (U) reduce_median();
        }

        if constexpr(ops == REDUCE_STD)
        {
            return (U) reduce_std();
        }

        if constexpr(ops == REDUCE_ARGMIN)
        {
            return (U) argmin();
        }

        if constexpr(ops == REDUCE_ARGMAX)
        {
            return (U) argmax();
        }
        panic("Invalid reduce operator\n");
    }

    T reduce_max() const {
        return reduce_ops<T, REDUCE_MAX>();
    }

    T reduce_min() const {
        return reduce_ops<T, REDUCE_MIN>();
    }

    template<typename U>
    U reduce_sum() const {
        return reduce_ops<U, REDUCE_SUM>();
    }

    double reduce_std() const 
    {
        if(_size == 0)
        {
            panic("Mean of empty array\n");
        }

        auto arr_square = square();
        double m2 = arr_square.reduce_mean();
        double m = reduce_mean();
        return ::sqrt(m2 - m*m);
    }

    T reduce_and() const {
        return reduce_ops<T, REDUCE_AND>();
    }

    T reduce_or() const {
        return reduce_ops<T, REDUCE_OR>();
    }

    T reduce_xor() const {
        return reduce_ops<T, REDUCE_XOR>();
    }

    int argmin() const {
        if(_size == 0)
        {
            panic("Min of empty array\n");
        }
        T min_val = _at(0);
        int imin = 0;
        for(int i = 1; i < _size; i++){
            T tmp = _at(i);
            if(tmp < min_val)
            {
                imin = i;
                min_val = tmp;
            }
        }
        return imin;
    }

    int argmax() const {
        if(_size == 0)
        {
            panic("Max of empty array\n");
        }
        T max_val = _at(0);
        int imax = 0;
        for(int i = 1; i < _size; i++){
            T tmp = _at(i);
            if(tmp > max_val)
            {
                imax = i;
                max_val = tmp;
            }
        }
        return imax;
    }

    void sort()
    {
        quickSort(get_ptr(), _step, NULL, 0, _size - 1);
    }

    Array<int> argsort() const
    {
        auto indexes = Array<int>::zeros(_size);
        for(int i=0; i < _size; i++) indexes[i] = i;
        auto tmp = copy();
        quickSort(tmp.get_ptr(), 1, indexes.get_ptr(), 0, _size-1);
        return indexes;
    }

    ~Array()
    {
        if(!_is_view && _data) {
            delete[] _data;
        }
    }
};

template<typename T>
static std::ostream & operator<<(std::ostream & os, const Array<T>& arr) 
{
    return print_list(os, arr);
}

template<typename INT>
Array<INT> slice_to_indexes(Slice idx, int size)
{
    idx = resolve_index(idx, size);
    assert(idx.step.value != 0);
    auto indexes = Array<INT>::zeros((idx.stop.value - idx.start.value)/idx.step.value);
    int n = 0;
    for(int i = idx.start.value; i != idx.stop.value; i += idx.step.value)
    {
        indexes[n] = i;
        n += 1;
    }
    return move(indexes);
}

template<typename BOOL, typename INT>
Array<INT> mask_to_indexes(const Array<BOOL>& masks)
{
    int n = 0;
    for(BOOL mask : masks) n += 1;
    auto indexes = Array<INT>::zeros(n);
    int idx = 0;

    for(int i = 0; i < masks.size(); i++)
    {
        if(masks[i])
        {
            indexes[idx] = i;
            idx += 1;
        }
    }

    return move(indexes);
}

template<typename T>
const T& _get_at(const T& val, int i)
{
    return val;
}

template<typename T>
const T& _get_at(const Array<T>& val, int i)
{
    return val._at(i);
}

#endif
