#include <vector>
#include <memory>
#include <cstddef>
#include <new>

template <class T, std::size_t Alignment>
class AlignedAllocator {
public:
    using value_type = T;

    AlignedAllocator() noexcept = default;

    template <class U>
    AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

    T* allocate(std::size_t n) {
        if (n == 0) return nullptr;

        void* ptr = ::operator new(
            n * sizeof(T),
            std::align_val_t{Alignment}
        );

        return static_cast<T*>(ptr);
    }

    void deallocate(T* p, std::size_t) noexcept {
        ::operator delete(
            p,
            std::align_val_t{Alignment}
        );
    }

    template <class U>
    struct rebind {
        using other = AlignedAllocator<U, Alignment>;
    };
};

template <class T1, std::size_t A1, class T2, std::size_t A2>
bool operator==(const AlignedAllocator<T1, A1>&,
                const AlignedAllocator<T2, A2>&) {
    return A1 == A2;
}

template <class T1, std::size_t A1, class T2, std::size_t A2>
bool operator!=(const AlignedAllocator<T1, A1>& a,
                const AlignedAllocator<T2, A2>& b) {
    return !(a == b);
}
