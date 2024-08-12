#ifndef __NP_H__
#define __NP_H__
#include "array4.h"

namespace np {
    template<typename T>
    Array<T> zeros(int n)
    {
        return Array<T>::zeros(n);
    }

    template<typename T>
    Array<T> ones(int n)
    {
        return Array<T>::ones(n);
    }

    template<typename T>
    Array<T> full(int n, T value)
    {
        return Array<T>::full(n, value);
    }

    template<typename T>
    Array<T> array(T* data, int n)
    {
        return Array<T>::new_array(data, n);
    }

    
    static Array<double> arange(double start, double stop, double step)
    {
        List<double> lst;
        if(step > 0)
        {
            while(start < stop)
            {
                lst.append(start);
                start += step;
            }
        }
        else if(step < 0) 
        {
            while(start > stop)
            {
                lst.append(start);
                start += step;
            }
        }
        else 
        {
            panic("Step 0 is now allowed");
        }
        return Array<double>::from_list(lst);
    }

    template<typename U, typename T>
    U sum(const Array<T>& arr)
    {
        return arr.template reduce_sum<U>();
    }

    template<typename T>
    T min(const Array<T>& arr)
    {
        return arr.reduce_min();
    }

    template<typename T>
    T max(const Array<T>& arr)
    {
        return arr.reduce_max();
    }

    template<typename T>
    double mean(const Array<T>& arr)
    {
        return arr.reduce_mean();
    }

    template<typename T>
    double median(const Array<T>& arr)
    {
        return arr.reduce_median();
    }

    template<typename T>
    double percentile(const Array<T>& arr, double p)
    {
        return arr.percentile(p);
    }

    template<typename T>
    double std_(const Array<T>& arr)
    {
        return arr.reduce_std();
    }

    template<typename BOOL>
    bool all(const Array<BOOL>& arr)
    {
        return arr.reduce_and();
    }

    template<typename BOOL>
    bool any_(const Array<BOOL>& arr)
    {
        return arr.reduce_or();
    }

    template<typename T>
    Array2<T> reshape(Array<T>& arr, Tuple<int, int> new_shape)
    {
        return arr.reshape(new_shape);
    }
 
    template<typename T>
    int argmin(const Array<T>& arr)
    {
        return arr.argmin();
    }

    template<typename T>
    int argmax(const Array<T>& arr)
    {
        return arr.argmax();
    }

    template<typename T>
    Array<T> sort(const Array<T>& arr)
    {
        auto res = arr.copy();
        res.sort();
        return res;
    }

    template<typename T>
    Array<int> argsort(const Array<T>& arr)
    {
        return arr.argsort();
    }

    template<typename T, typename Shape>
    auto zeros(Shape shape)
    {
        return full(shape, 0);
    }

    template<typename T, typename Shape>
    auto ones(Shape shape)
    {
        return full(shape, 1);
    }

    template<typename T, typename Shape>
    auto full(Shape shape, T value)
    {
        if constexpr(shape.dim() == 2)
            return Array2<T>::full(shape, value);

        else if constexpr(shape.dim() == 3)
            return Array3<T>::full(shape, value);
    }

    template<typename T, typename Shape>
    auto array(T data[], Shape shape)
    {
        if constexpr(shape.dim() == 2)
        {
            return Array2<T>::new_array2(data, shape); 
        }
        if constexpr(shape.dim() == 3)
        {
            return Array3<T>::new_array3(data, shape); 
        }
        if constexpr(shape.dim() == 4)
        {
            return Array4<T>::new_array4(data, shape); 
        }
    }

    template<typename T, typename ARR>
    Array<T> flatten(const ARR& arr)
    {
        Tuple<int, int> arr_shape = arr.shape();
        int n = arr.size();
        return Array<T>::from_raw(arr.copy().into_raw(), n);
    }

    template<typename T, typename ARR, typename Shape>
    auto reshape(const ARR& arr, Shape shape)
    {
        auto arr_shape = arr.shape();
        assert(shape.size() == arr_shape.size());
        return np::array(arr.copy().into_raw(), shape);
    }

    template<typename T>
    Array<T> min(const Array2<T>& arr, int ax)
    {
        return arr.reduce_min(ax);
    }

    template<typename T>
    Array<int> argmin(const Array2<T>&arr, int ax) 
    {
        return arr.argmin(ax); 
    }

    template<typename T>
    Array<T> max(const Array2<T>& arr, int ax)
    {
        return arr.reduce_max(ax);
    }

    template<typename T>
    Array<int> argmax(const Array2<T>&arr, int ax) 
    {
        return arr.argmax(ax); 
    }

    template<typename U, typename T>
    Array<U> sum(const Array2<T>& arr, int ax)
    {
        return arr.template reduce_sum<U>(ax);
    }

    template<typename T>
    Array<double> mean(const Array2<T>& arr, int ax)
    {
        return arr.reduce_mean(ax);
    }

    template<typename T>
    Array<double> median(const Array2<T>& arr, int ax)
    {
        return arr.reduce_median(ax);
    }

    template<typename T>
    Array<double> std_(const Array2<T>& arr, int ax)
    {
        return arr.reduce_std(ax);
    }

    template<typename T>
    Array<double> percentile(const Array2<T>& arr, double p, int ax)
    {
        return arr.percentile(p, ax);
    }

    template<typename T>
    Array<T> all(const Array2<T>& arr, int ax) 
    {
        return arr.reduce_and(ax);
    }

    template<typename T>
    Array<T> any_(const Array2<T>&arr, int ax)
    {
        return arr.reduce_or(ax);
    }

    template<typename T>
    T min(const Array2<T>& arr)
    {
        return arr.reduce_min();
    }

    template<typename T>
    int argmin(const Array2<T>& arr) 
    {
        return arr.argmin();
    }

    template<typename T>
    T max(const Array2<T>& arr)
    {
        return arr.reduce_max();
    }

    template<typename T>
    int argmax(const Array2<T>& arr) 
    {
        return arr.argmax();
    }

    template<typename U, typename T>
    U sum(const Array2<T>& arr)
    {
        return arr.template reduce_sum<U>();
    }

    template<typename T>
    double mean(const Array2<T>& arr)
    {
        return arr.reduce_mean();
    }

    template<typename T>
    double median(const Array2<T>& arr)
    {
        return arr.reduce_median();
    }

    template<typename T>
    double percentile(const Array2<T>& arr, double p) 
    {
        return arr.percentile(p);
    }

    template<typename T>
    double std_(const Array2<T>& arr)
    {
        return arr.reduce_std();
    }

    template<typename T>
    bool all(const Array2<T>& arr)
    {
        return arr.reduce_and();
    }

    template<typename T>
    bool any_(const Array2<T>& arr)
    {
        return arr.reduce_or();
    }

    template<typename T>
    Array2<T> sort(const Array2<T>& arr, int ax=1)
    {
        Array2<T> res = arr.copy();
        res.sort(ax);
        return res;
    }

    template<typename T>
    Array2<int> argsort(const Array2<T>& arr, int ax=1)
    {
        return arr.argsort(ax);
    }
}

#endif
