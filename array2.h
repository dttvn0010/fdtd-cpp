#ifndef __ARRAY2_H__
#define __ARRAY2_H__
#include "array.h"

template<typename T>
class Array2Iterator {
    int _step = 0;
    int _item_step = 0;
    int _item_length = 0;
    T* _data = NULL;
    const size_t* _p_signature = NULL;
public:
    Array2Iterator()
    {
        
    }
    Array2Iterator(T* data, int step, int item_step, int item_length, const size_t* p_signature){
        _step = step;
        _item_step = item_step;
        _item_length = item_length;
        _data = data;
        _p_signature = p_signature;
    }

    Array2Iterator& operator++(){
        _data += _step;
        return *this;
    }
    
    Array2Iterator& operator +=(int step){
        _data += _step * step;
        return *this;
    }

    bool operator != (Array2Iterator& rhs) {
        return _data != rhs._data;
    }

    Array<T> operator*() const {
        return Array<T>(
            _data,          //data
            _item_length,   // size
            _item_step,     // step
            true,           // is_view
            _p_signature    // _p_signature
        );
    }
};

struct CArray2 
{
    void* data;
    isize shape1;
    isize shape2;
    isize step1;
    isize step2;
    isize dtype;
};

template<typename T>
class Array2 : public Object
{
    template<typename U> friend class Array;
    template<typename U> friend class Array2;
    template<typename U> friend class Array3;
    template<typename U> friend class Array4;

protected:
    T* _data = NULL;
    int _step1 = 0;
    int _step2 = 0;

    int _shape1 = 0;
    int _shape2 = 0;
    bool _is_view = false;
    const size_t* _p_signature = NULL;

    T& _at(int i1, int i2) const
    {
        return _data[i2 * _step2 + i1 * _step1];
    }
    

    Array<T> _row_at(int i1) const
    {
        return Array<T>(
            _data + i1 * _step1,
            _shape2,
            _step2,
            true,
            get_p_signature()
        );
    }

    Array<T> _column_at(int i2) const
    {
        return Array<T>(
            _data + i2 * _step2,
            _shape1,
            _step1,
            true,
            get_p_signature()
        );
    }


    Array2(T* data, int shape1, int shape2, int step1, int step2, bool is_view=false, const size_t* p_signature=NULL)
    {
        _data = data;
        _shape1 = shape1;
        _shape2 = shape2;
        _step1 = step1;
        _step2 = step2;
        _is_view = is_view;
        _p_signature = p_signature;
        if(p_signature)
        {
            _signature = *_p_signature;
        }
    }

    void operator =(const Array2&) = delete;

    Array2(const Array2&) = delete;

    const size_t* get_p_signature() const
    {
        return _is_view? _p_signature : &_signature;
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
public:
    Array2()
    {

    }

    Array2(const CArray2* carr, char* p_signature): 
        Array2((T*) carr->data, carr->shape1, carr->shape2, carr->step1, carr->step2, true, p_signature)
    {

    }

    Array2(Array2&& oth)
    {
        oth._check_valid();
        _data = oth._data;
        _step1 = oth._step1;
        _step2 = oth._step2;
        _shape1 = oth._shape1;
        _shape2 = oth._shape2;
        _signature = oth._signature;
        _p_signature = oth._p_signature;
        _is_view = oth._is_view;

        oth._data = NULL;
        oth._shape1 = oth._shape2 = 0;
    }

    void operator =(Array2&& oth)
    {
        _check_valid();
        oth._check_valid();
        if(_data && !_is_view) delete[] _data;

        _data = oth._data;
        _shape1 = oth._shape1;
        _shape2 = oth._shape2;
        _step1 = oth._step1;
        _step2 = oth._step2;
        _is_view = oth._is_view;
        _signature = oth._signature;
        _p_signature = oth._p_signature;

        oth._data = NULL;
        oth._shape1 = 0;
        oth._shape2 = 0;
    }

    T& at(int i1_, int i2_) const
    {
        _check_valid();
        int i1 = i1_ < 0 ? _shape1 + i1_ : i1_;
        int i2 = i2_ < 0 ? _shape2 + i2_ : i2_;
        #ifdef DEBUG
        if(i1 < 0 || i1 >= _shape1)
        {
            panic("Index %d is out of range %d\n", i1_, _shape1);
        }
        if(i2 < 0 || i2 >= _shape2)
        {
            panic("Index %d is out of range %d\n", i2_, _shape2);
        }
        #endif
        return _at(i1, i2);
    }

    Array<T> get_row(int i1_) const
    {
        _check_valid();
        int i1 = i1_ < 0 ? _shape1 + i1_ : i1_;

        #ifdef DEBUG
        if(i1 < 0 || i1 >= _shape1)
        {
            panic("Index %d is out of range %d\n", i1_, _shape1);
        }
        #endif
        return _row_at(i1);
    }

    Array<T> get_column(int i2_) const
    {
        _check_valid();
        int i2 = i2_ < 0 ? _shape2 + i2_ : i2_;

        #ifdef DEBUG
        if(i2 < 0 || i2 >= _shape2)
        {
            panic("Index %d is out of range %d\n", i2_, _shape2);
        }
        #endif
        return _column_at(i2);
    }

    auto get(int i1, Slice i2)
    {
        return get_row(i1).get(i2);
    }

    static Array2 from_raw(T* raw, Tuple<int, int> shape)
    {
        int n = shape.size();
        return Array2(raw, shape.item1, shape.item2, shape.item2, 1);
    }

    static Array2 new_array2(T* raw, Tuple<int, int> shape)
    {
        int n = shape.size();
        T* data = new T[n];
        for(int i = 0; i < n; i++) data[i] = raw[i];
        return Array2(data, shape.item1, shape.item2, shape.item2, 1);
    }
    
    static Array2 full(Tuple<int, int> shape, T value)
    {
        int n = shape.item1 * shape.item2;
        T* data = new T[n];
        for(int i = 0; i < n; i++) data[i] = value;
        return Array2(data, shape.item1, shape.item2, shape.item2, 1);
    }

    static Array2 zeros(Tuple<int, int> shape)
    {
        return full(shape, (T) 0);
    }

    static Array2 ones(Tuple<int, int> shape)
    {
        return full(shape, (T) 1);
    }
    
    int size() const 
    {
        return _shape1;
    }

    Tuple<int, int> shape() const
    {
        _check_valid();
        return Tuple<int, int>(_shape1, _shape2);
    }

    Tuple<int, int> steps() const
    {
        _check_valid();
        return Tuple<int, int>(_step1, _step2);
    }

    bool is_view() const
    {
        return _is_view;
    }

    template<typename U>
    Array2<U> as_type() const {
        _check_valid();
        Array2<U> res = Array2<U>::zeros(shape());

        for(int i1 = 0; i1 < _shape1; i1++) 
        {
            res._row_at(i1).copy_from(_row_at(i1));
        }
        return res;
    }

    Array2<T> copy() const
    {
        return as_type<T>();
    }

    T* into_raw() 
    {
        _check_valid();
        T* data = _data;
        _data = NULL;
        return data;
    }

    template<typename INT>
    Array<INT> get_indexes(Slice idx, int ax)
    {
        int sz = ax == 0 ? _shape1 : _shape2;
        return slice_to_indexes<INT>(idx, sz); 
    }

    template<typename INT>
    Array<INT> get_indexes(const Array<INT>& indexes, int ax)
    {
        return indexes.copy();
        //return indexes.view();
    }

    template<typename BOOL, typename INT>
    Array<INT> get_indexes(const Array<BOOL>& masks, int ax)
    {
        return mask_to_indexes(masks);
    }

    template<typename Idx1, typename Idx2>
    void get_by_indexes_to(const Idx1& idx1, const Idx2& idx2, Array2<T>& out)
    {
        _check_valid();
        auto row_indexes = get_indexes(idx1, 0);
        auto col_indexes = get_indexes(idx2, 1);

        for(int i = 0; i < row_indexes.size(); i++)
        {
            auto out_row = out.get_row(i);
            get_row(row_indexes._at(i)).get_by_indexes_to(col_indexes, out_row);
        }
    }

    template<typename Idx1, typename Idx2> 
    Array2<T> get(const Idx1& idx1, const Idx2& idx2)
    {
        auto row_indexes = get_indexes(idx1, 0);
        auto col_indexes = get_by_indexes(idx2, 1);

        int n_row = row_indexes.size();
        int n_col = col_indexes.size();

        auto res = Array2<T>::zeros(Tuple<int, int>(n_row, n_col));
        get_by_indexes_to(row_indexes, col_indexes, res);
        return res;
    }

    template<typename Idx1>
    Array<T> get(const Idx1& idx1)
    {
        return get(idx1, Slice::all());
    }

    template<typename U, typename Idx1, typename Idx2>
    void set(const Idx1& idx1, const Idx2& idx2, const U& val)
    {
        _check_valid();
        auto row_indexes = get_indexes(idx1, 0);
        auto col_indexes = get_indexes(idx2, 1);

        for(int i = 0; i < row_indexes.size(); i++)
        {
            get_row(row_indexes._at(i)).set_by_indexes(col_indexes, val);
        }
    }

    template<typename Idx1, typename Idx2>
    void set(const Idx1& idx1, const Idx2& idx2, const Array2<T>& val)
    {
        _check_valid();
        auto row_indexes = get_indexes(idx1, 0);
        auto col_indexes = get_indexes(idx2, 1);

        for(int i = 0; i < row_indexes.size(); i++)
        {
            auto val_row = val.get_row(i);
            get_row(row_indexes._at(i)).set_by_indexes(col_indexes, val_row);
        }
    }

    template<typename U, typename Idx1>
    void set(const Idx1& idx1, const U& val){
        set(idx1, Slice::all(), val);
    }

    template<typename Idx1, typename Idx2>
    void get_by_pair_indexes_to(const Idx1& idx1, const Idx2& idx2, Array<T>& out)
    {
        _check_valid();
        auto row_indexes = get_indexes(idx1, 0);
        auto col_indexes = get_indexes(idx2, 1);

        int n_row = row_indexes.size();
        int n_col = col_indexes.size();
        int n_item = (n_row > n_col)? n_row : n_col;

        if(n_row > 1 && n_col > 1 && n_row != n_col)
        {
            panic("Row indexes size does not match column indexes size: %d vs %d\n", n_row, n_col);
        }

        for(int i = 0; i < n_item; i++)
        {
            int i1 = i < n_row? row_indexes._at(i) : row_indexes._at(0);
            int i2 = i < n_col? col_indexes._at(i) : col_indexes._at(0);
            out._at(i) = at(i1, i2);
        }
    }

    template<typename Idx1, typename Idx2>
    Array<T> get_by_pair_indexes(const Idx1& idx1, const Idx2& idx2)
    {
        auto row_indexes = get_indexes(idx1, 0);
        auto col_indexes = get_indexes(idx2, 1);

        int n_row = row_indexes.size();
        int n_col = col_indexes.size();
        int n_item = (n_row > n_col)? n_row : n_col;
        auto res = Array<T>::zeros(n_item);

        get_by_pair_indexes_to(row_indexes, col_indexes, res);
        return res;
    }

    template<typename Idx1, typename Idx2, typename U>
    void set_by_pair_indexes(const Idx1& idx1, const Idx2& idx2, const U& val)
    {
        auto row_indexes = get_indexes(idx1, 0);
        auto col_indexes = get_indexes(idx2, 1);

        int n_row = row_indexes.size();
        int n_col = col_indexes.size();

        int n_item = (n_row > n_col)? n_row : n_col;

        for(int i = 0; i < n_item; i++)
        {
            int i1 = i < n_row? row_indexes._at(i) : row_indexes._at(0);
            int i2 = i < n_col? col_indexes._at(i) : col_indexes._at(0);
            at(i1, i2) = _get_at(val, i);
        }
    }

    template<typename BOOL>
    void get_by_2d_masks_to(const Array2<BOOL>& masks, Array<T>& out)
    {
        _check_valid();
        assert(masks._shape1 == _shape1 && masks._shape2 == _shape2);
        auto indexes = mask_to_indexes<BOOL, int>(masks.flatten()); 
        flatten().get_by_indexes_to(indexes, out);
    }

    Array<T> get(const Array2<bool>& masks)
    {
        assert(masks._shape1 == _shape1 && masks._shape2 == _shape2);
        return flatten().get_by_masks(masks.flatten());
    }

    template<typename BOOL, typename U>
    int _set_by_masks(const Array<BOOL>& masks, int offset, const U& val)
    {
        int n = 0;
        for(auto mask: masks)
        {
            if(mask) n += 1;
        }

        //assert(n == val.size());

        int idx = 0;
        int i = 0;
        for(int i1 = 0; i1 < _shape1; i1++) 
        {
            auto row = get_row(i1);
            for(int i2 = 0; i2 < _shape2; i2++)
            {
                if(masks._at(i))
                {
                    row._at(i2) = _get_at(val, idx + offset);
                    idx += 1;
                }
                i += 1;
            }
        }
        return n;
    }

    template<typename BOOL, typename U>
    void set_by_2d_masks(const Array2<bool>& masks, const U& val)
    {
        _check_valid();
        assert(masks._shape1 == _shape1 && masks._shape2 == _shape2);
        set_by_2d_masks(masks.flatten(), 0, val);
    }

    Array<T> flatten() const
    {
        _check_valid();
        int n = _shape1 * _shape2;
        Array<T> res = Array<T>::zeros(n);
        int idx = 0;
        for(auto row : *this)
        {
            for(int i2 = 0; i2 < _shape2; i2++)
            {
                res._at(idx) = row._at(i2);
                idx ++; 
            }
        }
        return res;
    }

    void transpose_to(Array2<T>& arr) {
        _check_valid();
        arr._check_valid();

        for(int i2=0; i2 < _shape2; i2++){
            arr._row_at(i2).copy_from(_column_at(i2));
        }
    }

    Array2<T> transpose() {
        _check_valid();
        auto arr = Array2::zeros(Tuple<int, int>(_shape2, _shape1));
        transpose_to(arr);
        return arr;
    }

    void set_all(const T& value)
    {
        _check_valid();
        for(int i1=0; i1 < _shape1; i1++){
            _row_at(i1).set_all(value);
        }
    }

    void assign(const T& value)
    {
        set_all(value);
    }

    template<typename U>
    void set_all_rows(const Array<U>& row) {
        _check_valid();
        row._check_valid();

        if(row._size != _shape2)
        {
            panic("Array sizes mismatch : %d vs %d\n", _shape2, row._size);
        }
        for(int i1 = 0; i1 < _shape1; i1++)
        {
            _row_at(i1).copy_from(row);
        }
    }

    template<typename U>
    void assign(const Array<U>& row)
    {
        set_all_rows(row);
    }

    template<typename U>
    void copy_from(const Array2<U>& oth)
    {
        _check_valid();
        oth._check_valid();

        if(_shape1 != oth._shape1 || _shape2 != oth._shape2)
        {
            panic("Array size mismatch (%d, %d) vs (%d, %d)\n", _shape1, _shape2, oth._shape1, oth._shape2);
        }

        for(int i1 = 0; i1 < _shape1; i1++) 
        {
            _row_at(i1).copy_from(oth._row_at(i1));
        }
    }

    template<typename U>
    void assign(const Array2<U>& oth)
    {
        copy_from(oth);
    }

    template<typename U>
    void set_all_columns(const Array<U>& col) 
    {
        _check_valid();
        col._check_valid();

        if(col._size != _shape1)
        {
            panic("Array sizes mismatch: %d vs %d\n", _shape1, col._size);
        }
        for(int i2 = 0; i2 < _shape2; i2++)
        {
            _column_at(i2).copy_from(col);
        }
    }

    Array2Iterator<T> begin() const
    {
        _check_valid();
        return Array2Iterator<T>(_data, _step1, _step2, _shape2, get_p_signature());
    }

    Array2Iterator<T> end() const
    {
        _check_valid();
        return Array2Iterator<T>(_data + _step1 * _shape1, _step1, _step2, _shape2, get_p_signature());
    }
    
    Array<T> operator[](int idx)
    {
        return get_row(idx);
    }

    Array2<T> get(Slice idx1, Slice idx2)
    {
        _check_valid();
        idx1 = resolve_index(idx1, _shape1);
        idx2 = resolve_index(idx2, _shape2);

        int start1 = idx1.start.value;
        int stop1 = idx1.stop.value;
        int step1 = idx1.step.value;

        int start2 = idx2.start.value;
        int stop2 = idx2.stop.value;
        int step2 = idx2.step.value;

        assert(step1 != 0 && step2 != 0);

        return Array2<T>(
            _data + _step1 * start1 + _step2 * start2,      // _data
            (stop1 - start1)/step1,                         // _shape1
            (stop2 - start2)/step2,                         // _shape2
            _step1 * step1,                                 // _step1
            _step2 * step2,                                 // _step2
            true,
            get_p_signature()
        );
    }

    Array2<T> get(Slice idx1)
    {
        _check_valid();
        idx1 = resolve_index(idx1, _shape1);

        int start1 = idx1.start.value;
        int stop1 = idx1.stop.value;
        int step1 = idx1.step.value;

        assert(step1 != 0);

        return Array2<T>(
            _data + _step1 * start1,                        // _data
            (stop1 - start1)/step1,                         // _shape1
            _shape2,                                        // _shape2
            _step1 * step1,                                 // _step1
            _step2,                                         // _step2
            true,
            get_p_signature()
        );
    }

    Array<T> get(Slice idx1, int idx2)
    {
        _check_valid();

        idx1 = resolve_index(idx1, _shape1);
        int start1 = idx1.start.value;
        int stop1 = idx1.stop.value;
        int step1 = idx1.step.value;

        assert(step1 != 0);

        return Array<T>(
            _data + _step1 * start1 + _step2 * idx2,        // _data
            (stop1 - start1)/step1,                         // _size
            _step1 * step1,                                 // _step
            true,
            get_p_signature()
        );
    }

    template<typename Idx1, typename U>
    void set(Idx1 idx1, const U& val)
    {
        get(idx1).assign(val);
    }

    template<typename Idx1, typename Idx2, typename U>
    void set(Idx1 idx1, Idx2 idx2, const U& val)
    {
        get(idx1, idx2).assign(val);
    }

    template<int func>
    void unary_ops_compute(Array2<T>& res) const
    {
        _check_valid();
        res._check_valid();
        assert(func <= NUM_FUNC);

        for(int i1 = 0; i1 < _shape1; i1++)
        {
            auto res_row = res._row_at(i1); 
            _row_at(i1).template unary_ops_compute<func>(res_row);
        }
    }

    template<int func>
    Array2<T> unary_ops() const
    {
        auto res = Array2<T>::zeros(shape()) ; 
        unary_ops_compute<func>(res); 
        return res;
    }

    Array2<T> operator ~() {
        return unary_ops<FUNC_NOT>();
    }

    Array2<T> operator -() {
        return unary_ops<FUNC_NEG>();
    }

    Array2<T> abs()  const {
        return unary_ops<FUNC_ABS>();
    }

    Array2<T> round() const {
        return unary_ops<FUNC_ROUND>();
    }

    Array2<T> ceil() const {
        return unary_ops<FUNC_CEIL>();
    }

    Array2<T> floor() const {
        return unary_ops<FUNC_FLOOR>();
    }

    Array2<T> square() const {
        return unary_ops<FUNC_SQUARE>();
    }

    Array2<T> sqrt() const {
        return unary_ops<FUNC_SQRT>();
    }

    Array2<T> exp() const {
        return unary_ops<FUNC_EXP>();
    }

    Array2<T> log() const {
        return unary_ops<FUNC_LOG>();
    }

    Array2<T> sin() const {
        return unary_ops<FUNC_SIN>();
    }

    Array2<T> cos() const {
        return unary_ops<FUNC_COS>();
    }

    Array2<T> tan() const {
        return unary_ops<FUNC_TAN>();
    }

    Array2<T> asin() const {
        return unary_ops<FUNC_ASIN>();
    }

    Array2<T> acos() const {
        return unary_ops<FUNC_ACOS>();
    }

    Array2<T> atan() const {
        return unary_ops<FUNC_ATAN>();
    }

    template<int ops, typename U>
    void bin_ops_assign(const U& rhs) {
        _check_valid();

        if constexpr(ops > NUM_OP)
        {
            panic("Invalid operator: %d\n", ops);
        }

        for(auto row : *this)
        {
            row.template bin_ops_assign<ops>(rhs);
        }
    }

    template<int ops>
    void bin_ops_assign(const Array2<T>& rhs) {
        _check_valid();
        rhs._check_valid();

        if constexpr(ops > NUM_OP)
        {
            panic("Invalid operator: %d\n", ops);
        }

        if(
            (rhs._shape1 != _shape1 && rhs._shape1 != 1) ||
            (rhs._shape2 != _shape2 && rhs._shape2 != 1)
        )
        {
            panic("Dimensions mismatch: (%d,%d) vs (%d, %d)", _shape1, _shape2, rhs._shape1, rhs._shape2);
        }

        if(rhs._shape1 == 1)
        {
            for(auto row : *this)
            {
                row.template bin_ops_assign<ops>(rhs._row_at(0));
            }
        }
        else
        {
            for(auto [l_row, r_row] : zip(*this, rhs))
            {
                l_row.template bin_ops_assign<ops>(r_row);    
            }
        }
    }

    template<typename U>
    void operator +=(const U& rhs) {
        return bin_ops_assign<OP_ADD>(rhs);
    }

    template<typename U>
    void operator -=(const U& rhs) {
        return bin_ops_assign<OP_SUB>(rhs);
    }

    template<typename U>
    void operator *=(const U& rhs) {
        return bin_ops_assign<OP_MUL>(rhs);
    }

    template<typename U>
    void operator /=(const U& rhs) {
        return bin_ops_assign<OP_DIV>(rhs);
    }

    template<typename U>
    void operator %=(const U& rhs) {
        return bin_ops_assign<OP_MOD>(rhs);
    }

    template<typename U>
    void operator &=(const U& rhs) {
        bin_ops_assign<OP_AND>(rhs);
    }

    template<typename U> 
    void operator |=(const U& rhs) {
        bin_ops_assign<OP_OR>(rhs);
    }

    template<typename U>
    void operator ^=(const U& rhs) {
        bin_ops_assign<OP_XOR>(rhs);
    }

    template<typename U, int ops, typename V>
    static void bin_ops_compute(const Array2<T>& lhs, const V& rhs, Array2<U>& res)
    {
        lhs._check_valid();

        if constexpr(ops > NUM_OP) 
        {
            panic("Invalid operator: %d", ops);
        }

        for(auto[l_row, res_row] : zip(lhs, res))
        {
            Array<T>::template bin_ops_compute<U, ops>(
                l_row,
                rhs,
                res_row
            );
        }
    }

    template<typename U, int ops, typename V>
    Array2<U> bin_ops(const V& rhs) const {
        auto res = Array2<U>::zeros(shape());
        bin_ops_compute<U, ops>(*this, rhs, res); 
        return res;
    }

    template<typename U, int ops>
    static void bin_ops_compute(const Array2<T>& lhs, 
            const Array2<T>& rhs, 
            Array2<U>& res)
    {
        lhs._check_valid();
        rhs._check_valid();

        if constexpr(ops > NUM_OP) 
        {
            panic("Invalid operator: %d", ops);
        }

        int shape1 = lhs._shape1 > rhs._shape1 ? lhs._shape1 : rhs._shape1;

        if(lhs._shape1 > 1 && rhs._shape1 > 1)
        {
            for(auto [l_row,r_row,res_row] : zip(lhs, rhs, res))
            {
                Array<T>::template bin_ops_compute<U, ops>(
                    l_row,
                    r_row,
                    res_row
                );
            }
        }
        else
        {
            for(int i1 = 0; i1 < shape1; i1++)
            {
                auto l_row = lhs._row_at(lhs._shape1 > 1? i1 : 0);
                auto r_row = rhs._row_at(rhs._shape1 > 1? i1:0);
                auto res_row = res._row_at(i1);

                Array<T>::template bin_ops_compute<U, ops>(
                    l_row,
                    r_row,
                    res_row
                );
            }
        }
    }

    template<typename U, int ops>
    Array2<U> bin_ops(const Array2<T>& rhs) const
    {
        if(
            (rhs._shape1 != _shape1 && rhs._shape1 != 1 && _shape1 != 1) ||
            (rhs._shape2 != _shape2 && rhs._shape2 != 1 && _shape2 != 1)
        )
        {
            panic("Dimensions mismatch: (%d,%d) vs (%d, %d)", _shape1, _shape2, rhs._shape1, rhs._shape2);
        }

        int shape1 = _shape1 > rhs._shape1 ? _shape1 : rhs._shape1;
        int shape2 = _shape2 > rhs._shape2 ? _shape2 : rhs._shape2;

        auto res = Array2<U>::zeros(Tuple<int, int>(shape1, shape2));
        bin_ops_compute<U, ops>(*this, rhs, res);
        return res;
    }

    template<typename V>
    Array2<T> operator +(const V& rhs) const {
        return bin_ops<T, OP_ADD>(rhs);
    }

    template<typename V>
    Array2<T> operator -(const V& rhs) const{
        return bin_ops<T, OP_SUB>(rhs);
    }

    template<typename V>
    Array2<T> operator *(const V& rhs) const{
        return bin_ops<T, OP_MUL>(rhs);
    }

    template<typename V>
    Array2<T> operator /(const V& rhs) const{
        return bin_ops<T, OP_DIV>(rhs);
    }

    template<typename V>
    Array2<T> operator %(const V& rhs) const{
        return bin_ops<T, OP_MOD>(rhs);
    }

    template<typename V>
    Array2<T> pow(const V& rhs) const {
        return bin_ops<T, OP_POW>(rhs);
    }

    template<typename V>
    Array2<T> operator &(const V& rhs) const{
        return bin_ops<T, OP_AND>(rhs);
    }

    template<typename V>
    Array2<T> operator |(const V& rhs) const{
        return bin_ops<T, OP_OR>(rhs);
    }

    template<typename V>
    Array2<T> operator ^(const V& rhs) const{
        return bin_ops<T, OP_XOR>(rhs);
    }

    template<typename V>
    Array2<bool> operator >(const V& rhs) const{
        return bin_ops<bool, OP_GT>(rhs);
    }

    template<typename V>
    Array2<bool> operator >=(const V& rhs) const{
        return bin_ops<bool, OP_GTE>(rhs);
    }

    template<typename V>
    Array2<bool> operator <(const V& rhs) const{
        return bin_ops<bool, OP_LT>(rhs);
    }

    template<typename V>
    Array2<bool> operator <=(const V& rhs) const{
        return bin_ops<bool, OP_LTE>(rhs);
    }

    template<typename V>
    Array2<bool> operator ==(const V& rhs) const{
        return bin_ops<bool, OP_EQ>(rhs);
    }

    template<typename V>
    Array2<bool> operator !=(const V& rhs) const{
        return bin_ops<bool, V, OP_NE>(rhs);
    }

    template<typename U, int ops>
    void reduce_ops_compute(int ax, Array<U>& res) const 
    {
        _check_valid();

        if(ax > 1)
        {
            panic("Invalid axis: %d\n", ax);
        }

        int sz = ax == 0? _shape2 : _shape1;

        if(ax == 0)
        {
            for(int i = 0; i < sz; i++) 
                res._at(i) = get_column(i).template reduce_ops<U, ops>();
        }else
        {
            for(int i = 0; i < sz; i++)
                res._at(i) = get_row(i).template reduce_ops<U, ops>();
        }

    }

    template<typename U, int ops>
    Array<U> reduce_ops(int ax) const {
        int sz = ax == 0? _shape2 : _shape1;
        auto res = Array<U>::zeros(sz); 
        reduce_ops_compute<U, ops>(ax, res);
        return res;
    }

    Array<double> percentile(double p, int ax) const {
        _check_valid();

        if(ax > 1)
        {
            panic("Invalid axis: %d\n", ax);
        }

        int sz = ax == 0? _shape2 : _shape1;

        auto res = Array<double>::zeros(sz);
        if(ax == 0)
        {
            for(int i = 0; i < sz; i++) 
                res._at(i) = get_column(i).percentile(p);
        }else
        {
            for(int i = 0; i < sz; i++)
                res._at(i) = get_row(i).percentile(p);
        }

        return res;
    }

    Array<T> reduce_max(int ax) const {
        return reduce_ops<T, REDUCE_MAX>(ax);
    }

    Array<T> reduce_min(int ax) const {
        return reduce_ops<T, REDUCE_MIN>(ax);
    }

    template<typename U>
    Array<U> reduce_sum(int ax) const {
        return reduce_ops<U, REDUCE_SUM>(ax);
    }

    Array<double> reduce_mean(int ax) const {
        return reduce_ops<double, REDUCE_MEAN>(ax);
    }

    Array<double> reduce_median(int ax) const {
        return reduce_ops<double, REDUCE_MEDIAN>(ax);
    }

    Array<double> reduce_std(int ax) const {
        return reduce_ops<double, REDUCE_STD>(ax);
    }

    Array<T> reduce_and(int ax) const {
        return reduce_ops<T, REDUCE_AND>(ax);
    }

    Array<T> reduce_or(int ax) const {
        return reduce_ops<T, REDUCE_OR>(ax);
    }

    Array<T> reduce_xor(int ax) const {
        return reduce_ops<T, REDUCE_XOR>(ax);
    }

    Array<int> argmax(int ax) const
    {
        return reduce_ops<int, REDUCE_ARGMAX>(ax);
    }

    Array<int> argmin(int ax) const
    {
        return reduce_ops<int, REDUCE_ARGMIN>(ax);
    }

    T reduce_max() const{
        return flatten().reduce_max();
    }

    T reduce_min() const{
        return flatten().reduce_min();
    }

    
    T reduce_sum() const {
        return flatten().template reduce_sum<T>();
    }

    double reduce_mean() const {
        return flatten().reduce_mean();
    }

    double reduce_median() const {
        return flatten().reduce_median();
    }

    double percentile(double p) {
        return flatten().percentile(p);
    }

    double reduce_std() const {
        return flatten().reduce_std();
    }

    T reduce_and() const {
        return flatten().reduce_and();
    }

    T reduce_or() const {
        return flatten().reduce_or();
    }

    T reduce_xor() const {
        return flatten().reduce_xor();
    }

    void sort(int ax=1)
    {
        _check_valid();
        if(ax > 1)
        {
            panic("Invalid axis: %d\n", ax);
        }

        if(ax == 1)
        {
            for(int i = 0; i < _shape1; i++)
            {
                _row_at(i).sort();
            }
        } else
        {
            for(int i = 0; i < _shape2; i++)
            {
                _column_at(i).sort();
            }
        }
    }

    Array2<int> argsort(int ax=1) const
    {
        _check_valid();
        if(ax > 1)
        {
            panic("Invalid axis: %d\n", ax);
        }
        auto arr = Array2<int>::zeros(shape());
        if(ax == 1) 
        {
            for(int i = 0; i < _shape1; i++)
            {
                arr._row_at(i).copy_from(_row_at(i).argsort());
            }
        } else
        {
            for(int i = 0; i < _shape2; i++)
            {
                arr._column_at(i).copy_from(_column_at(i).argsort());
            }
        }
        return arr;
    }

    int argmax() const
    {
        return flatten().argmax();
    }

    int argmin() const
    {
        return flatten().argmin();
    }

    ~Array2()
    {
        if(!_is_view && _data) {
            delete[] _data;
        }
    }
};

template<typename T>
static std::ostream & operator<<(std::ostream & os, const Array2<T>& arr) 
{
    return print_list(os, arr);
}

#endif
