#ifndef __ARRAY4_H__
#define __ARRAY4_H__
#include "array.h"
#include "array3.h"
#include <cassert>

template<typename T>
class Array4 : public Object
{
    template<typename U> friend class Array;
    template<typename U> friend class Array2;
    template<typename U> friend class Array3;
    template<typename U> friend class Array4;
protected:
    T* _data = NULL;
    int _step1 = 0;
    int _step2 = 0;
    int _step3 = 0;
    int _step4 = 0;

    int _shape1 = 0;
    int _shape2 = 0;
    int _shape3 = 0;
    int _shape4 = 0;

    bool _is_view = false;
    const size_t* _p_signature = NULL;

    T& _at(int i1, int i2, int i3, int i4) const
    {
        return _data[i4 * _step4 + i3 * _step3 + i2 * _step2 + i1 * _step1];
    }

    Array4(T* data, int shape1, int shape2, int shape3, int shape4, 
        int step1, int step2, int step3, int step4,
        bool is_view=false, const size_t* p_signature=NULL)
    {
        _data = data;
        _shape1 = shape1;
        _shape2 = shape2;
        _shape3 = shape3;
        _shape4 = shape4;
        _step1 = step1;
        _step2 = step2;
        _step3 = step3;
        _step4 = step4;
        _is_view = is_view;
        _p_signature = p_signature;
        if(p_signature)
        {
            _signature = *_p_signature;
        }
    }

    void operator =(const Array4&) = delete;

    Array4(const Array4&) = delete;

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

    Array3<T> _get(int idx) const
    {
        return Array3<T>(
            _data + idx * _step1,
            _shape2,
            _shape3,
            _shape4,
            _step2,
            _step3,
            _step4,
            true,
            get_p_signature()
        );
    }

public:
    Array4()
    {

    }
    Array4(Array4&& oth)
    {
        oth._check_valid();
        _data = oth._data;
        _step1 = oth._step1;
        _step2 = oth._step2;
        _step3 = oth._step3;
        _step4 = oth._step4;
        _shape1 = oth._shape1;
        _shape2 = oth._shape2;
        _shape3 = oth._shape3;
        _shape4 = oth._shape4;
        _signature = oth._signature;
        _p_signature = oth._p_signature;
        _is_view = oth._is_view;

        oth._data = NULL;
        oth._shape1 = oth._shape2 = oth._shape3 = oth._shape4 = 0;
    }

    void operator =(Array4&& oth)
    {
        _check_valid();
        oth._check_valid();
        if(_data && !_is_view) delete[] _data;

        _data = oth._data;
        _shape1 = oth._shape1;
        _shape2 = oth._shape2;
        _shape3 = oth._shape3;
        _shape4 = oth._shape4;
        _step1 = oth._step1;
        _step2 = oth._step2;
        _step3 = oth._step3;
        _step4 = oth._step4;
        _is_view = oth._is_view;
        _signature = oth._signature;
        _p_signature = oth._p_signature;

        oth._data = NULL;
        oth._shape1 = 0;
        oth._shape2 = 0;
        oth._shape3 = 0;
        oth._shape4 = 0;
    }

    static Array4 new_array4(T* raw, Tuple<int, int, int, int> shape)
    {
        int n = shape.size();
        T* data = new T[n];
        for(int i = 0; i < n; i++) data[i] = raw[i];
        return Array4(
            data, 
            shape.item1,                                    // shape1
            shape.item2,                                    // shape2
            shape.item3,                                    // shape3
            shape.item4,                                    // shape4
            shape.item4 * shape.item3 * shape.item2,        // step1
            shape.item4 * shape.item3,                      // step2
            shape.item4,                                    // step3
            1                                               // step4
        );
    }

    static Array4 from_raw(T* raw, Tuple<int, int, int, int> shape)
    {
        int n = shape.size();
        return Array4(
            raw, 
            shape.item1,                                    // shape1
            shape.item2,                                    // shape2
            shape.item3,                                    // shape3
            shape.item4,                                    // shape4
            shape.item4 * shape.item3 * shape.item2,        // step1
            shape.item4 * shape.item3,                      // step2
            shape.item4,                                    // step3
            1                                               // step4
        );
    }

    static Array4 full(Tuple<int, int, int, int> shape, T value)
    {
        int n = shape.size();
        T* data = new T[n];
        for(int i = 0; i < n; i++) data[i] = value;
        return Array4(
            data, 
            shape.item1,                                    // shape1
            shape.item2,                                    // shape2
            shape.item3,                                    // shape3
            shape.item4,                                    // shape4
            shape.item4 * shape.item3 * shape.item2,        // step1
            shape.item4 * shape.item3,                      // step2
            shape.item4,                                    // step3
            1                                               // step4
        );
    }

    static Array4 zeros(Tuple<int, int, int, int> shape)
    {
        return full(shape, (T) 0);
    }

    static Array4 ones(Tuple<int, int, int, int> shape)
    {
        return full(shape, (T) 1);
    }

    Tuple<int, int, int, int> shape() const
    {
        _check_valid();
        return Tuple<int, int, int, int>(_shape1, _shape2, _shape3, _shape4);
    }

    T& at(int i1_, int i2_, int i3_, int i4_) const
    {
        _check_valid();
        int i1 = i1_ < 0 ? _shape1 + i1_ : i1_;
        int i2 = i2_ < 0 ? _shape2 + i2_ : i2_;
        int i3 = i3_ < 0 ? _shape3 + i3_ : i3_;
        int i4 = i4_ < 0 ? _shape4 + i4_ : i4_;

        #ifdef DEBUG
        if(i1 < 0 || i1 >= _shape1)
        {
            panic("Index %d is out of range %d\n", i1_, _shape1);
        }
        if(i2 < 0 || i2 >= _shape2)
        {
            panic("Index %d is out of range %d\n", i2_, _shape2);
        }
        if(i3 < 0 || i3 >= _shape3)
        {
            panic("Index %d is out of range %d\n", i3_, _shape3);
        }
        if(i4 < 0 || i4 >= _shape4)
        {
            panic("Index %d is out of range %d\n", i4_, _shape4);
        }
        #endif
        return _at(i1, i2, i3, i4);
    }

    Array2<T> get(const List<int>& idx1, const List<int>& idx2, const List<int>& idx3) const
    {
        assert(idx1.size() == idx2.size() && idx2.size() == idx3.size());
        auto result = Array2<T>::zeros(Tuple<int, int>(idx1.size(), _shape4));
        for(int i = 0; i < (int) idx1.size(); i++)
        {
            for(int j = 0; j < _shape4; j++)
            {
                result.at(i,j) = _at(idx1[i], idx2[i], idx3[i], j);
            }
        }
        return result;
    }

    Array2<T> get(int i1_, Slice idx2, Slice idx3, int i4_) const
    {
        _check_valid();
        int i1 = i1_ < 0 ? _shape1 + i1_ : i1_;
        idx2 = resolve_index(idx2, _shape2);
        idx3 = resolve_index(idx3, _shape3);
        int i4 = i4_ < 0 ? _shape4 + i4_ : i4_;

        int start2 = idx2.start.value;
        int stop2 = idx2.stop.value;
        int step2 = idx2.step.value;

        int start3 = idx3.start.value;
        int stop3 = idx3.stop.value;
        int step3 = idx3.step.value;

        assert(step2 != 0 && step3 != 0);

        #ifdef DEBUG
        if(i1 < 0 || i1 >= _shape1)
        {
            panic("Index %d is out of range %d\n", i1_, _shape1);
        }
        if(i4 < 0 || i4 >= _shape4)
        {
            panic("Index %d is out of range %d\n", i4_, _shape4);
        }
        #endif

        return Array2<T>(
            _data + _step1 * i1 + _step2 * start2 + _step3 * start3 + _step4 * i4,      // _data
            (stop2 - start2)/step2,                                         // _shape2
            (stop3 - start3)/step3,                                         // _shape3
            _step2 * step2,                                                 // _step2
            _step3 * step3,                                                 // _step3
            true,
            get_p_signature()
        );
    }

    Array2<T> get(Slice idx1, int i2_, Slice idx3, int i4_) const
    {
        _check_valid();
        idx1 = resolve_index(idx1, _shape1);
        idx3 = resolve_index(idx3, _shape3);

        int i2 = i2_ < 0 ? _shape2 + i2_ : i2_;
        int i4 = i4_ < 0 ? _shape4 + i4_ : i4_;
        
        int start1 = idx1.start.value;
        int stop1 = idx1.stop.value;
        int step1 = idx1.step.value;

        int start3 = idx3.start.value;
        int stop3 = idx3.stop.value;
        int step3 = idx3.step.value;

        assert(step1 != 0 && step3 != 0);

        #ifdef DEBUG
        if(i2 < 0 || i2 >= _shape2)
        {
            panic("Index %d is out of range %d\n", i2_, _shape2);
        }
        if(i4 < 0 || i4 >= _shape4)
        {
            panic("Index %d is out of range %d\n", i4_, _shape4);
        }
        #endif

        return Array2<T>(
            _data + _step1 * start1 + _step2 * i2 + _step3 * start3 + _step4 * i4_,      // _data
            (stop1 - start1)/step1,                                         // _shape1
            (stop3 - start3)/step3,                                         // _shape3
            _step1 * step1,                                                 // _step1
            _step3 * step3,                                                 // _step3
            true,
            get_p_signature()
        );
    }

    Array2<T> get(Slice idx1, Slice idx2, int i3_, int i4_) const
    {
        _check_valid();
        idx1 = resolve_index(idx1, _shape1);
        idx2 = resolve_index(idx2, _shape2);
        int i3 = i3_ < 0 ? _shape3 + i3_ : i3_;
        int i4 = i4_ < 0 ? _shape4 + i4_ : i4_;

        int start1 = idx1.start.value;
        int stop1 = idx1.stop.value;
        int step1 = idx1.step.value;

        int start2 = idx2.start.value;
        int stop2 = idx2.stop.value;
        int step2 = idx2.step.value;


        assert(step1 != 0 && step2 != 0);

        #ifdef DEBUG
        if(i3 < 0 || i3 >= _shape3)
        {
            panic("Index %d is out of range %d\n", i3_, _shape3);
        }
        if(i4 < 0 || i4 >= _shape4)
        {
            panic("Index %d is out of range %d\n", i4_, _shape4);
        }
        #endif

        return Array2<T>(
            _data + _step1 * start1 + _step2 * start2 + _step3 * i3 + _step4 * i4,      // _data
            (stop1 - start1)/step1,                                         // _shape1
            (stop2 - start2)/step2,                                         // _shape2
            _step1 * step1,                                                 // _step1
            _step2 * step2,                                                 // _step2
            true,
            get_p_signature()
        );
    }

    Array3<T> operator[](int idx) const
    {
        _check_valid();
        int i1 = idx < 0 ? _shape1 + idx : idx;

        #ifdef DEBUG
        if(i1 < 0 || i1 >= _shape1)
        {
            panic("Index %d is out of range %d\n", idx, _shape1);
        }
        #endif
        return Array3<T>(
            _data + i1 * _step1,
            _shape2,
            _shape3,
            _shape4,
            _step2,
            _step3,
            _step4,
            true,
            get_p_signature()
        );
    }

    Array3<T> get(Slice idx1, int i2_) const
    {
        _check_valid();

        int i2 = i2_ < 0 ? _shape4 + i2_ : i2_;
        #ifdef DEBUG
        if(i2 < 0 || i2 >= _shape2)
        {
            panic("Index %d is out of range %d\n", i2_, _shape2);
        }
        #endif

        idx1 = resolve_index(idx1, _shape1);
        int start1 = idx1.start.value;
        int stop1 = idx1.stop.value;
        int step1 = idx1.step.value;
    
        assert(step1 != 0);

        return Array3<T>(
            _data + _step1 * start1 + _step2 * i2,                          // _data
            (stop1 - start1)/step1,                                         // _shape1
            _shape3,                                                        // _shape2
            _shape4,                                                        // _shape3
            _step1 * step1,                                                 // _step1
            _step3,                                                         // _step2
            _step4,                                                         // _step3
            true,
            get_p_signature()
        );
    }

    Array3<T> get(Slice idx1, Slice idx2, int i3_) const
    {
        _check_valid();

        int i3 = i3_ < 0 ? _shape3 + i3_ : i3_;
        #ifdef DEBUG
        if(i3 < 0 || i3 >= _shape3)
        {
            panic("Index %d is out of range %d\n", i3_, _shape3);
        }
        #endif

        idx1 = resolve_index(idx1, _shape1);
        idx2 = resolve_index(idx2, _shape2);

        int start1 = idx1.start.value;
        int stop1 = idx1.stop.value;
        int step1 = idx1.step.value;

        int start2 = idx2.start.value;
        int stop2 = idx2.stop.value;
        int step2 = idx2.step.value;

        assert(step1 != 0 && step2 != 0);

        return Array3<T>(
            _data + _step1 * start1 + _step2 * start2 + _step3 * i3,        // _data
            (stop1 - start1)/step1,                                         // _shape1
            (stop2 - start2)/step2,                                         // _shape2
            _shape4,                                                        // _shape3
            _step1 * step1,                                                 // _step1
            _step2 * step2,                                                 // _step2
            _step4,                                                         // _step3
            true,
            get_p_signature()
        );
    }

    Array4<T> get(Slice idx1, Slice idx2, Slice idx3, Slice idx4) const
    {
        _check_valid();
        idx1 = resolve_index(idx1, _shape1);
        idx2 = resolve_index(idx2, _shape2);
        idx3 = resolve_index(idx3, _shape3);
        idx4 = resolve_index(idx4, _shape4);

        int start1 = idx1.start.value;
        int stop1 = idx1.stop.value;
        int step1 = idx1.step.value;

        int start2 = idx2.start.value;
        int stop2 = idx2.stop.value;
        int step2 = idx2.step.value;

        int start3 = idx3.start.value;
        int stop3 = idx3.stop.value;
        int step3 = idx3.step.value;

        int start4 = idx4.start.value;
        int stop4 = idx4.stop.value;
        int step4 = idx4.step.value;

        assert(step1 != 0 && step2 != 0 && step3 != 0 && step4 != 0);

        return Array4<T>(
            _data + _step1 * start1 + _step2 * start2 + _step3 * start3 + _step4 * start4,      // _data
            (stop1 - start1)/step1,                                         // _shape1
            (stop2 - start2)/step2,                                         // _shape2
            (stop3 - start3)/step3,                                         // _shape3
            (stop4 - start4)/step4,                                         // _shape4
            _step1 * step1,                                                 // _step1
            _step2 * step2,                                                 // _step2
            _step3 * step3,                                                 // _step3
            _step4 * step4,                                                 // _step3
            true,
            get_p_signature()
        );
    }

    Array4<T> get(Tuple<Slice, Slice, Slice, Slice> idx) const
    {
        return get(idx.item1, idx.item2, idx.item3, idx.item4);
    }

    Array3<T> get(Slice idx1, Slice idx2, Slice idx3, int i4_) const
    {
        _check_valid();

        int i4 = i4_ < 0 ? _shape4 + i4_ : i4_;
        #ifdef DEBUG
        if(i4 < 0 || i4 >= _shape4)
        {
            panic("Index %d is out of range %d\n", i4_, _shape4);
        }
        #endif

        idx1 = resolve_index(idx1, _shape1);
        idx2 = resolve_index(idx2, _shape2);
        idx3 = resolve_index(idx3, _shape3);

        int start1 = idx1.start.value;
        int stop1 = idx1.stop.value;
        int step1 = idx1.step.value;

        int start2 = idx2.start.value;
        int stop2 = idx2.stop.value;
        int step2 = idx2.step.value;

        int start3 = idx3.start.value;
        int stop3 = idx3.stop.value;
        int step3 = idx3.step.value;

        assert(step1 != 0 && step2 != 0 && step3 != 0);

        return Array3<T>(
            _data + _step1 * start1 + _step2 * start2 + _step3 * start3 + _step4 * i4,      // _data
            (stop1 - start1)/step1,                                         // _shape1
            (stop2 - start2)/step2,                                         // _shape2
            (stop3 - start3)/step3,                                         // _shape3
            _step1 * step1,                                                 // _step1
            _step2 * step2,                                                 // _step2
            _step3 * step3,                                                 // _step3
            true,
            get_p_signature()
        );
    }

    Array3<T> get(Tuple<Slice, Slice, Slice, int> idx) const
    {
        return get(idx.item1, idx.item2, idx.item3, idx.item4);
    }

    template<typename U>
    Array4<U> as_type() const {
        _check_valid();
        auto res = Array4<U>::zeros(shape());

        for(int i1 = 0; i1 < _shape1; i1++) 
        {
            res._get(i1).copy_from(_get(i1));
        }
        return res;
    }

    Array4<T> copy() const
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

    template<typename U>
    void copy_from(const Array4<U>& oth)
    {
        _check_valid();
        oth._check_valid();

        if(_shape1 != oth._shape1 || _shape2 != oth._shape2 || _shape3 != oth._shape3 || _shape4 != oth._shape4)
        {
            panic("Array size mismatch (%d, %d, %d, %d) vs (%d, %d, %d, %d)\n", 
                _shape1, _shape2, _shape3, _shape4,
                oth._shape1, oth._shape2, oth._shape3, oth._shape4);
        }

        for(int i1 = 0; i1 < _shape1; i1++) 
        {
            _get(i1).copy_from(oth._get(i1));
        }
    }

    template<typename U>
    void assign(const Array4<U>& oth)
    {
        copy_from(oth);
    }

    void set_all(const T& value)
    {
        _check_valid();
        for(int i1=0; i1 < _shape1; i1++){
            _get(i1).set_all(value);
        }
    }

    template<int func>
    void unary_ops_compute(Array4<T>& res) const
    {
        _check_valid();
        res._check_valid();
        assert(func <= NUM_FUNC);

        for(int i1 = 0; i1 < _shape1; i1++)
        {
            auto res_row = res._get(i1); 
            _get(i1).template unary_ops_compute<func>(res_row);
        }
    }

    template<int func>
    Array4<T> unary_ops() const
    {
        auto res = Array4<T>::zeros(shape()) ; 
        unary_ops_compute<func>(res); 
        return res;
    }

    Array4<T> operator ~() {
        return unary_ops<FUNC_NOT>();
    }

    Array4<T> operator -() {
        return unary_ops<FUNC_NEG>();
    }

    Array4<T> abs()  const {
        return unary_ops<FUNC_ABS>();
    }

    Array4<T> round() const {
        return unary_ops<FUNC_ROUND>();
    }

    Array4<T> ceil() const {
        return unary_ops<FUNC_CEIL>();
    }

    Array4<T> floor() const {
        return unary_ops<FUNC_FLOOR>();
    }

    Array4<T> square() const {
        return unary_ops<FUNC_SQUARE>();
    }

    Array4<T> sqrt() const {
        return unary_ops<FUNC_SQRT>();
    }

    Array4<T> exp() const {
        return unary_ops<FUNC_EXP>();
    }

    Array4<T> log() const {
        return unary_ops<FUNC_LOG>();
    }

    Array4<T> sin() const {
        return unary_ops<FUNC_SIN>();
    }

    Array4<T> cos() const {
        return unary_ops<FUNC_COS>();
    }

    Array4<T> tan() const {
        return unary_ops<FUNC_TAN>();
    }

    Array4<T> asin() const {
        return unary_ops<FUNC_ASIN>();
    }

    Array4<T> acos() const {
        return unary_ops<FUNC_ACOS>();
    }

    Array4<T> atan() const {
        return unary_ops<FUNC_ATAN>();
    }

    template<int ops, typename U>
    void bin_ops_assign(const U& rhs) {
        _check_valid();

        if constexpr(ops > NUM_OP)
        {
            panic("Invalid operator: %d\n", ops);
        }

        for(int i = 0; i < _shape1; i++)
        {
            _get(i).template bin_ops_assign<ops>(rhs);
        }
    }

    template<int ops>
    void bin_ops_assign(const Array4<T>& rhs) {
        _check_valid();
        rhs._check_valid();

        if constexpr(ops > NUM_OP)
        {
            panic("Invalid operator: %d\n", ops);
        }

        if(
            (rhs._shape1 != _shape1 && rhs._shape1 != 1) ||
            (rhs._shape2 != _shape2 && rhs._shape2 != 1) ||
            (rhs._shape3 != _shape3 && rhs._shape3 != 1) ||
            (rhs._shape4 != _shape4 && rhs._shape4 != 1)
        )
        {
            panic("Dimensions mismatch: (%d, %d, %d, %d) vs (%d, %d, %d, %d)", 
                _shape1, _shape2, _shape3, _shape4,
                rhs._shape1, rhs._shape2, rhs._shape3, rhs._shape4
            );
        }

        if(rhs._shape1 == 1)
        {
            for(int i1 = 0; i1 < _shape1; i1++)
            {
                _get(i1).template bin_ops_assign<ops>(rhs._get(0));
            }
        }
        else
        {
            for(int i1 = 0; i1 < _shape1; i1++)
            {
                _get(i1).template bin_ops_assign<ops>(rhs._get(i1));
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


    template<typename U, int ops, typename V>
    static void bin_ops_compute(const Array4<T>& lhs, const V& rhs, Array4<U>& res)
    {
        lhs._check_valid();

        if constexpr(ops > NUM_OP) 
        {
            panic("Invalid operator: %d", ops);
        }

        for(int i1 = 0; i1 < lhs._shape1; i1++)
        {
            auto tmp = res._get(i1);
            Array3<T>::template bin_ops_compute<U, ops, V>(
                lhs._get(i1),
                rhs,
                tmp
            );
        }
    }

    template<typename U, int ops, typename V>
    Array4<U> bin_ops(const V& rhs) const {
        auto res = Array4<U>::zeros(shape());
        bin_ops_compute<U, ops>(*this, rhs, res); 
        return res;
    }

    template<typename U, int ops>
    static void bin_ops_compute(const Array4<T>& lhs, 
            const Array4<T>& rhs, 
            Array4<U>& res)
    {
        lhs._check_valid();
        rhs._check_valid();

        if constexpr(ops > NUM_OP) 
        {
            panic("Invalid operator: %d", ops);
        }

        int shape1 = lhs._shape1 > rhs._shape1 ? lhs._shape1 : rhs._shape1;
        assert(lhs._shape1 == 1 || lhs._shape1 == shape1);
        assert(rhs._shape1 == 1 || rhs._shape1 == shape1);

        if(lhs._shape1 > 1 && rhs._shape1 > 1)
        {
            for(int i1 = 0; i1 < shape1; i1++)
            {
                auto tmp = res._get(i1);
                Array3<T>::template bin_ops_compute<U, ops>(
                    lhs._get(i1),
                    rhs._get(i1),
                    tmp
                );
            }            
        }
        else
        {
            for(int i1 = 0; i1 < shape1; i1++)
            {
                auto tmp = res._get(i1);
                Array3<T>::template bin_ops_compute<U, ops>(
                    lhs._get(lhs._shape1 > 1? i1 : 0),
                    rhs._get(rhs._shape1 > 1? i1:0),
                    tmp
                );
            }
        }
    }

    template<typename U, int ops>
    Array4<U> bin_ops(const Array4<T>& rhs) const
    {
        if(
            (rhs._shape1 != _shape1 && rhs._shape1 != 1 && _shape1 != 1) ||
            (rhs._shape2 != _shape2 && rhs._shape2 != 1 && _shape2 != 1) ||
            (rhs._shape3 != _shape3 && rhs._shape3 != 1 && _shape3 != 1) ||
            (rhs._shape4 != _shape4 && rhs._shape4 != 1 && _shape4 != 1)
        )
        {
            panic("Dimensions mismatch: (%d, %d, %d, %d) vs (%d, %d, %d, %d)", 
                _shape1, _shape2, _shape3, _shape4,
                rhs._shape1, rhs._shape2, rhs._shape3, rhs._shape4
            );
        }

        int shape1 = _shape1 > rhs._shape1 ? _shape1 : rhs._shape1;
        int shape2 = _shape2 > rhs._shape2 ? _shape2 : rhs._shape2;
        int shape3 = _shape3 > rhs._shape3 ? _shape3 : rhs._shape3;
        int shape4 = _shape4 > rhs._shape4 ? _shape4 : rhs._shape4;

        auto res = Array4<U>::zeros(Tuple<int, int, int, int>(shape1, shape2, shape3 ,_shape4));
        bin_ops_compute<U, ops>(*this, rhs, res);
        return res;
    }

    template<typename V>
    Array4<T> operator +(const V& rhs) const {
        return bin_ops<T, OP_ADD>(rhs);
    }

    template<typename V>
    Array4<T> operator -(const V& rhs) const{
        return bin_ops<T, OP_SUB>(rhs);
    }

    template<typename V>
    Array4<T> operator *(const V& rhs) const{
        return bin_ops<T, OP_MUL>(rhs);
    }

    template<typename V>
    Array4<T> operator /(const V& rhs) const{
        return bin_ops<T, OP_DIV>(rhs);
    }

    ~Array4()
    {
        if(!_is_view && _data) {
            delete[] _data;
        }
    }
};
#endif