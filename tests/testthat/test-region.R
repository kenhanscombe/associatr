context("region")

test_that(
    "region() stops with a msg if no access token supplied",
    {
        expect_error(
            region(),
            "^Argument `token` requires an LDlinkR 'Personal Access Token'.*"
        )
    }
)
