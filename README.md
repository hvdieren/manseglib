# mantissa-segmentation

Repository for the MEng project "A runtime library for mantissa segmentation".

## Notes about current implementation
- can perform basic mantissa segmentation with expected outcomes

When compiled with -O4:
- performance is less than desired:
    - using heads only is slightly slower than double performance
    - using heads+tails (pairs) is about twice as slow as double


# ToDo List:
* Investigate and fix performance issues
    * too many conversions?
    * inefficient arithmetic operations?
    * ...

* Find larger and more "real-world" test data for PageRank
* Explore other algorithms