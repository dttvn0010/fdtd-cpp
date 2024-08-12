#ifndef __ARRAY3_H__
#define __ARRAY3_H__
#include "array2.h"

template<typename T>
class Array3 : public Object
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

    int _shape1 = 0;
    int _shape2 = 0;
    int _shape3 = 0;
    bool _is_view = false;
    const size_t* _p_signature = NULL;

    T& _at(int i1, int i2, int i3) const
    {
        return _data[i3 * _step3 + i2 * _step2 + i1 * _step1];
    }

    Array3(T* data, int shape1, int shape2, int shape3, int step1, int step2, int step3, 
        bool is_view=false, const size_t* p_signature=NULL)
    {
        _data = data;
        _shape1 = shape1;
        _shape2 = shape2;
        _shape3 = shape3;
        _step1 = step1;
        _step2 = step2;
        _step3 = step3;
        _is_view = is_view;
        _p_signature = p_signature;
        if(p_signature)
        {
            _signature = *_p_signature;
        }
    }

    void operator =(const Array3&) = delete;

    Array3(const Array3&) = delete;

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

    template<int ax>
    Array2<T> _get(int idx) const
    {
        if constexpr(ax==0)
        {
            return Array2<T>(
                _data + idx * _step1,
                _shape2,
                _shape3,
                _step2,
                _step3,
                true,
                get_p_signature()
            );
        }
        if constexpr(ax==1)
        {
            return Array2<T>(
                _data + idx * _step2,
                _shape1,
                _shape3,
                _step1,
                _step3,
                true,
                get_p_signature()
            );
        }
        if constexpr(ax==2)
        {
            return Array2<T>(
                _data + idx * _step3,
                _shape1,
                _shape2,
                _step1,
                _step2,
                true,
                get_p_signature()
            );
        }
    }

    int get_min_axis() const
    {
        if(_shape1 < _shape2)
        {
            return _shape1 < _shape3 ? 0 : 2;
        }
        else 
        {
            return _shape2 < _shape3 ? 1 : 2;
        }
    }

    bool is_compact() const
    {
        return _step3 == 1 && _step2 == _shape3 && _step1 == _shape2 * _shape3;
    }

    Array<T> get_compact_view() const
    {
        return Array<T>(_data, _shape1 * _shape2 * _shape3, 1, true, get_p_signature());
    }
public:
    Array3()
    {

    }

    Array3(Array3&& oth)
    {
        oth._check_valid();
        _data = oth._data;
        _step1 = oth._step1;
        _step2 = oth._step2;
        _step3 = oth._step3;
        _shape1 = oth._shape1;
        _shape2 = oth._shape2;
        _shape3 = oth._shape3;
        _signature = oth._signature;
        _p_signature = oth._p_signature;
        _is_view = oth._is_view;

        oth._data = NULL;
        oth._shape1 = oth._shape2 = oth._shape3 = 0;
    }

    void operator =(Array3&& oth)
    {
        _check_valid();
        oth._check_valid();
        if(_data && !_is_view) delete[] _data;

        _data = oth._data;
        _shape1 = oth._shape1;
        _shape2 = oth._shape2;
        _shape3 = oth._shape3;
        _step1 = oth._step1;
        _step2 = oth._step2;
        _step3 = oth._step3;
        _is_view = oth._is_view;
        _signature = oth._signature;
        _p_signature = oth._p_signature;

        oth._data = NULL;
        oth._shape1 = 0;
        oth._shape2 = 0;
        oth._shape3 = 0;
    }

    static Array3 new_array3(T* raw, Tuple<int, int, int> shape)
    {
        int n = shape.size();
        T* data = new T[n];
        for(int i = 0; i < n; i++) data[i] = raw[i];
        return Array3(
            data, 
            shape.item1,                        // shape1
            shape.item2,                        // shape2
            shape.item3,                        // shape3
            shape.item3 * shape.item2,          // step1
            shape.item3,                        // step2
            1                                   // step3
        );
    }

    static Array3 from_raw(T* raw, Tuple<int, int, int> shape)
    {
        int n = shape.size();
        return Array3(
            raw, 
            shape.item1,                        // shape1
            shape.item2,                        // shape2
            shape.item3,                        // shape3
            shape.item3 * shape.item2,          // step1
            shape.item3,                        // step2
            1                                   // step3
        );
    }
    
    static Array3 full(Tuple<int, int, int> shape, T value)
    {
        int n = shape.size();
        T* data = new T[n];
        for(int i = 0; i < n; i++) data[i] = value;
        return Array3(
            data, 
            shape.item1,                        // shape1
            shape.item2,                        // shape2
            shape.item3,                        // shape3
            shape.item3 * shape.item2,          // step1
            shape.item3,                        // step2
            1                                   // step3
        );
    }

    static Array3 zeros(Tuple<int, int, int> shape)
    {
        return full(shape, (T) 0);
    }

    static Array3 ones(Tuple<int, int, int> shape)
    {
        return full(shape, (T) 1);
    }

    Tuple<int, int, int> shape() const
    {
        _check_valid();
        return Tuple<int, int, int>(_shape1, _shape2, _shape3);
    }
   
    T& at(int i1_, int i2_, int i3_) const
    {
        _check_valid();
        int i1 = i1_ < 0 ? _shape1 + i1_ : i1_;
        int i2 = i2_ < 0 ? _shape2 + i2_ : i2_;
        int i3 = i3_ < 0 ? _shape3 + i3_ : i3_;

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
        #endif
        return _at(i1, i2, i3);
    }

    Array2<T> operator[](int idx) const
    {
        _check_valid();
        int i1 = idx < 0 ? _shape1 + idx : idx;

        #ifdef DEBUG
        if(i1 < 0 || i1 >= _shape1)
        {
            panic("Index %d is out of range %d\n", idx, _shape1);
        }
        #endif

        return _get<0>(i1);
    }

    Array3<T> get(Slice idx1, Slice idx2, Slice idx3) const
    {
        _check_valid();
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
            _data + _step1 * start1 + _step2 * start2 + _step3 * start3,      // _data
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

    Array<T> get(Slice idx1, int idx2, int idx3) const
    {
        _check_valid();
        idx1 = resolve_index(idx1, _shape1);
        int i2 = idx2 < 0 ? _shape2 + idx2 : idx2;
        int i3 = idx3 < 0 ? _shape3 + idx3 : idx3;

        #ifdef DEBUG
        if(i2 < 0 || i2 >= _shape2)
        {
            panic("Index %d is out of range %d\n", idx2, _shape2);
        }
        if(i3 < 0 || i3 >= _shape3)
        {
            panic("Index %d is out of range %d\n", idx3, _shape3);
        }
        #endif

        int start1 = idx1.start.value;
        int stop1 = idx1.stop.value;
        int step1 = idx1.step.value;

        assert(step1 != 0);

        return Array<T>(
            _data + _step1 * start1 + _step2 * i2 + _step3 * i3,            // _data
            (stop1 - start1)/step1,                                         // _shape1
            _step1 * step1,                                                 // _step1
            true,
            get_p_signature()
        );
    }

    Array<T> get(int idx1, Slice idx2, int idx3) const
    {
        _check_valid();
        idx2 = resolve_index(idx2, _shape2);
        int i1 = idx1 < 0 ? _shape1 + idx1 : idx1;
        int i3 = idx3 < 0 ? _shape3 + idx3 : idx3;

        #ifdef DEBUG
        if(i1 < 0 || i1 >= _shape1)
        {
            panic("Index %d is out of range %d\n", idx1, _shape1);
        }
        if(i3 < 0 || i3 >= _shape3)
        {
            panic("Index %d is out of range %d\n", idx3, _shape3);
        }
        #endif

        int start2 = idx2.start.value;
        int stop2 = idx2.stop.value;
        int step2 = idx2.step.value;

        assert(step2 != 0);

        return Array<T>(
            _data + _step1 * i1 + _step2 * start2 + _step3 * i3,            // _data
            (stop2 - start2)/step2,                                         // _shape1
            _step2 * step2,                                                 // _step1
            true,
            get_p_signature()
        );
    }

    Array<T> get(int idx1, int idx2, Slice idx3) const
    {
        _check_valid();
        idx3 = resolve_index(idx3, _shape3);
        int i1 = idx1 < 0 ? _shape1 + idx1 : idx1;
        int i2 = idx2 < 0 ? _shape2 + idx2 : idx2;

        #ifdef DEBUG
        if(i1 < 0 || i1 >= _shape1)
        {
            panic("Index %d is out of range %d\n", idx1, _shape1);
        }
        if(i2 < 0 || i2 >= _shape2)
        {
            panic("Index %d is out of range %d\n", idx2, _shape2);
        }
        #endif

        int start3 = idx3.start.value;
        int stop3 = idx3.stop.value;
        int step3 = idx3.step.value;

        assert(step3 != 0);

        return Array<T>(
            _data + _step1 * i1 + _step2 * i2 + _step3 * start3,            // _data
            (stop3 - start3)/step3,                                         // _shape1
            _step3 * step3,                                                 // _step1
            true,
            get_p_signature()
        );
    }

    template<typename U>
    Array3<U> as_type() const {
        _check_valid();
        auto res = Array3<U>::zeros(shape());
        res.copy_from(*this);
        return res;
    }

    Array3<T> copy() const
    {
        return as_type<T>();
    }

    T* into_raw() 
    {
        _check_valid();
        T* data = _data;
        _shape1 = _shape2 = _shape3 = 0;
        _data = NULL;
        return data;
    }

    template<typename U>
    void copy_from(const Array3<U>& oth)
    {
        _check_valid();
        oth._check_valid();

        if(_shape1 != oth._shape1 || _shape2 != oth._shape2 || _shape3 != oth._shape3)
        {
            panic("Array size mismatch (%d, %d, %d) vs (%d, %d, %d)\n", 
                _shape1, _shape2, _shape3,
                oth._shape1, oth._shape2, oth._shape3
            );
        }

        if(is_compact() && oth.is_compact())
        {
            get_compact_view().copy_from(oth.get_compact_view());
            return;
        }

        int min_ax = get_min_axis();
        if(min_ax == 0)
        {
            for(int i = 0; i < _shape1; i++) 
            {
                _get<0>(i).copy_from(oth.template _get<0>(i));
            }
        }
        if(min_ax == 1)
        {
            for(int i = 0; i < _shape2; i++) 
            {
                _get<1>(i).copy_from(oth.template _get<1>(i));
            }
        }
        if(min_ax == 2)
        {
            for(int i = 0; i < _shape3; i++) 
            {
                _get<2>(i).copy_from(oth.template _get<2>(i));
            }
        }
    }

    template<typename U>
    void assign(const Array3<U>& oth)
    {
        copy_from(oth);
    }

    void set_all(const T& value)
    {
        _check_valid();
        if(is_compact())
        {
            get_compact_view().set_all(value);
            return;
        }

        int min_ax = get_min_axis();
        if(min_ax == 0)
        {
            for(int i=0; i < _shape1; i++)
            {
                _get<0>(i).set_all(value);
            }
        }
        if(min_ax == 1)
        {
            for(int i=0; i < _shape2; i++)
            {
                _get<1>(i).set_all(value);
            }
        }
        if(min_ax == 2)
        {
            for(int i=0; i < _shape3; i++)
            {
                _get<2>(i).set_all(value);
            }
        }
    }

    template<int func>
    void unary_ops_compute(Array3<T>& res) const
    {
        _check_valid();
        res._check_valid();
        assert(func <= NUM_FUNC);

        if(is_compact() && res.is_compact())
        {
            auto tmp = res.get_compact_view();
            get_compact_view().template unary_ops_compute<func>(tmp);
        }
        else
        {
            int min_ax = get_min_axis();
            if(min_ax == 0)
            {
                for(int i1 = 0; i1 < _shape1; i1++)
                {
                    auto tmp = res.template _get<0>(i1); 
                    _get<0>(i1).template unary_ops_compute<func>(tmp);
                }
            }
            else if(min_ax == 1)
            {
                for(int i2 = 0; i2 < _shape2; i2++)
                {
                    auto tmp = res.template _get<1>(i2); 
                    _get<1>(i2).template unary_ops_compute<func>(tmp);
                }
            }
            else if(min_ax == 2)
            {
                for(int i3 = 0; i3 < _shape3; i3++)
                {
                    auto tmp = res.template _get<2>(i3); 
                    _get<2>(i3).template unary_ops_compute<func>(tmp);
                }
            }
        }
    }

    template<int func>
    Array3<T> unary_ops() const
    {
        auto res = Array3<T>::zeros(shape()) ; 
        unary_ops_compute<func>(res); 
        return res;
    }

    Array3<T> operator ~() {
        return unary_ops<FUNC_NOT>();
    }

    Array3<T> operator -() {
        return unary_ops<FUNC_NEG>();
    }

    Array3<T> abs()  const {
        return unary_ops<FUNC_ABS>();
    }

    Array3<T> round() const {
        return unary_ops<FUNC_ROUND>();
    }

    Array3<T> ceil() const {
        return unary_ops<FUNC_CEIL>();
    }

    Array3<T> floor() const {
        return unary_ops<FUNC_FLOOR>();
    }

    Array3<T> square() const {
        return unary_ops<FUNC_SQUARE>();
    }

    Array3<T> sqrt() const {
        return unary_ops<FUNC_SQRT>();
    }

    Array3<T> exp() const {
        return unary_ops<FUNC_EXP>();
    }

    Array3<T> log() const {
        return unary_ops<FUNC_LOG>();
    }

    Array3<T> sin() const {
        return unary_ops<FUNC_SIN>();
    }

    Array3<T> cos() const {
        return unary_ops<FUNC_COS>();
    }

    Array3<T> tan() const {
        return unary_ops<FUNC_TAN>();
    }

    Array3<T> asin() const {
        return unary_ops<FUNC_ASIN>();
    }

    Array3<T> acos() const {
        return unary_ops<FUNC_ACOS>();
    }

    Array3<T> atan() const {
        return unary_ops<FUNC_ATAN>();
    }

    template<int ops>
    void bin_ops_assign(const T& rhs) {
        _check_valid();

        if constexpr(ops > NUM_OP)
        {
            panic("Invalid operator: %d\n", ops);
        }

        if(is_compact())
        {
            get_compact_view().template bin_ops_assign<ops>(rhs);
            return;
        }

        int min_ax = get_min_axis();
        if(min_ax == 0)
        {
            for(int i1 = 0; i1 < _shape1; i1++)
            {
                _get<0>(i1).template bin_ops_assign<ops>(rhs);
            }
        }
        else if(min_ax == 1)
        {
            for(int i2 = 0; i2 < _shape2; i2++)
            {
                _get<1>(i2).template bin_ops_assign<ops>(rhs);
            }
        }
        else if(min_ax == 2)
        {
            for(int i3 = 0; i3 < _shape3; i3++)
            {
                _get<2>(i3).template bin_ops_assign<ops>(rhs);
            }
        }
    }

    template<int ops>
    void bin_ops_assign(const Array3<T>& rhs) {
        _check_valid();
        rhs._check_valid();

        if constexpr(ops > NUM_OP)
        {
            panic("Invalid operator: %d\n", ops);
        }

        if(
            (rhs._shape1 != _shape1 && rhs._shape1 != 1) ||
            (rhs._shape2 != _shape2 && rhs._shape2 != 1) ||
            (rhs._shape3 != _shape3 && rhs._shape3 != 1)
        )
        {
            panic("Dimensions mismatch: (%d, %d, %d) vs (%d, %d, %d)", 
                _shape1, _shape2, _shape3,
                rhs._shape1, rhs._shape2, rhs._shape3
            );
        }

        if(is_compact() && rhs.is_compact() && _shape1 == rhs._shape1 && _shape2 == rhs._shape2 && _shape3 == rhs._shape3)
        {
            get_compact_view().template bin_ops_assign<ops>(rhs.get_compact_view());
        }
        else
        {
            int min_ax = get_min_axis();
            if(min_ax == 0)
            {
                for(int i1 = 0; i1 < _shape1; i1++)
                {
                    _get<0>(i1).template bin_ops_assign<ops>(rhs.template _get<0>(i1 < rhs._shape1 ? i1 : 0));
                }
            }
            else if(min_ax == 1)
            {
                for(int i2 = 0; i2 < _shape2; i2++)
                {
                    _get<1>(i2).template bin_ops_assign<ops>(rhs.template _get<1>(i2 < rhs._shape2 ? i2 : 0));
                }
            }
            else if(min_ax == 2)
            {
                for(int i3 = 0; i3 < _shape3; i3++)
                {
                    _get<2>(i3).template bin_ops_assign<ops>(rhs.template _get<2>(i3 < rhs._shape3 ? i3 : 0));
                }
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

    template<typename U, int ops>
    static void bin_ops_compute(const Array3<T>& lhs, const T& rhs, Array3<U>& res)
    {
        lhs._check_valid();

        if constexpr(ops > NUM_OP) 
        {
            panic("Invalid operator: %d", ops);
        }

        if(lhs.is_compact() && res.is_compact())
        {
            auto tmp = res.get_compact_view();
            Array<T>::template bin_ops_compute<U, ops>(lhs.get_compact_view(), rhs, tmp);
            return;
        }

        int min_ax = lhs.get_min_axis();

        if(min_ax == 0)
        {
            for(int i = 0; i < lhs._shape1; i++)
            {
                auto tmp = res.template _get<0>(i);
                Array2<T>::template bin_ops_compute<U, ops>(
                    lhs.template _get<0>(i),
                    rhs,
                    tmp
                );
            }
        }

        if(min_ax == 1)
        {
            for(int i = 0; i < lhs._shape2; i++)
            {
                auto tmp = res.template _get<1>(i);
                Array2<T>::template bin_ops_compute<U, ops>(
                    lhs.template _get<1>(i),
                    rhs,
                    tmp
                );
            }
        }

        if(min_ax == 2)
        {
            for(int i = 0; i < lhs._shape3; i++)
            {
                auto tmp = res.template _get<2>(i);
                Array2<T>::template bin_ops_compute<U, ops>(
                    lhs.template _get<2>(i),
                    rhs,
                    tmp
                );
            }
        }
    }

    template<typename U, int ops, typename V>
    Array3<U> bin_ops(const V& rhs) const {
        auto res = Array3<U>::zeros(shape());
        bin_ops_compute<U, ops>(*this, rhs, res); 
        return res;
    }

    template<typename U, int ops>
    static void bin_ops_compute(const Array3<T>& lhs, 
            const Array3<T>& rhs, 
            Array3<U>& res)
    {
        lhs._check_valid();
        rhs._check_valid();

        if constexpr(ops > NUM_OP) 
        {
            panic("Invalid operator: %d", ops);
        }

        int shape1 = lhs._shape1 > rhs._shape1 ? lhs._shape1 : rhs._shape1;
        int shape2 = lhs._shape2 > rhs._shape2 ? lhs._shape2 : rhs._shape2;
        int shape3 = lhs._shape3 > rhs._shape3 ? lhs._shape3 : rhs._shape3;

        if(lhs.is_compact() && rhs.is_compact() && res.is_compact() && 
            lhs._shape1 == rhs._shape1 && lhs._shape2 == rhs._shape2 && lhs._shape3 == rhs._shape3)
        {
            auto tmp = res.get_compact_view();
            Array<T>::template bin_ops_compute<U, ops>(
                lhs.get_compact_view(),
                rhs.get_compact_view(),
                tmp
            );
        }
        else
        {
            int min_ax = lhs.get_min_axis();
            if(min_ax == 0)
            {
                for(int i1 = 0; i1 < shape1; i1++)
                {
                    auto tmp = res.template _get<0>(i1);

                    Array2<T>::template bin_ops_compute<U, ops>(
                        lhs.template _get<0>(i1 < lhs._shape1 ? i1 : 0),
                        rhs.template _get<0>(i1 < rhs._shape1 ? i1 : 0),
                        tmp
                    );
                }
            }
            else if(min_ax == 1)
            {
                for(int i2 = 0; i2 < shape2; i2++)
                {
                    auto tmp = res.template _get<1>(i2);

                    Array2<T>::template bin_ops_compute<U, ops>(
                        lhs.template _get<1>(i2 < lhs._shape2 ? i2 : 0),
                        rhs.template _get<1>(i2 < rhs._shape2 ? i2 : 0),
                        tmp
                    );
                }
            }
            else if(min_ax == 2)
            {
                for(int i3 = 0; i3 < shape3; i3++)
                {
                    auto tmp = res.template _get<2>(i3);

                    Array2<T>::template bin_ops_compute<U, ops>(
                        lhs.template _get<2>(i3 < lhs._shape3 ? i3 : 0),
                        rhs.template _get<2>(i3 < rhs._shape3 ? i3 : 0),
                        tmp
                    );
                }
            }
        }
    }

    template<typename U, int ops>
    Array3<U> bin_ops(const Array3<T>& rhs) const
    {
        if(
            (rhs._shape1 != _shape1 && rhs._shape1 != 1 && _shape1 != 1) ||
            (rhs._shape2 != _shape2 && rhs._shape2 != 1 && _shape2 != 1) ||
            (rhs._shape3 != _shape3 && rhs._shape3 != 1 && _shape3 != 1)
        )
        {
            panic("Dimensions mismatch: (%d, %d, %d) vs (%d, %d, %d)", 
                _shape1, _shape2, _shape3,
                rhs._shape1, rhs._shape2, rhs._shape3
            );
        }

        int shape1 = _shape1 > rhs._shape1 ? _shape1 : rhs._shape1;
        int shape2 = _shape2 > rhs._shape2 ? _shape2 : rhs._shape2;
        int shape3 = _shape3 > rhs._shape3 ? _shape3 : rhs._shape3;

        auto res = Array3<U>::zeros(Tuple<int, int, int>(shape1, shape2, shape3));
        bin_ops_compute<U, ops>(*this, rhs, res);
        return res;
    }

    template<typename V>
    Array3<T> operator +(const V& rhs) const {
        return bin_ops<T, OP_ADD>(rhs);
    }

    template<typename V>
    Array3<T> operator -(const V& rhs) const{
        return bin_ops<T, OP_SUB>(rhs);
    }

    template<typename V>
    Array3<T> operator *(const V& rhs) const{
        return bin_ops<T, OP_MUL>(rhs);
    }

    template<typename V>
    Array3<T> operator /(const V& rhs) const{
        return bin_ops<T, OP_DIV>(rhs);
    }

    ~Array3()
    {
        if(!_is_view && _data) {
            delete[] _data;
        }
    }
};
#endif