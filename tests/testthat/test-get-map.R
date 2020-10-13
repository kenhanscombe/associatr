context("get_map")

test_that(
    "get_map stops and prints error for invalid pop",
    {
        expect_error(
            get_map("NON"),
            "^Invalid population code.*"
        )
    }
)
