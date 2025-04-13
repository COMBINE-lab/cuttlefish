cc_library(
    name = "epoch",
    hdrs = ["epoch.h"],
    deps = [
        "@parlaylib//parlay:primitives",
    ],
)

cc_library(
    name = "lock",
    hdrs = ["lock.h"],
)

cc_library(
    name = "unordered_map",
    hdrs = ["unordered_map.h"],
    deps = [
        ":epoch",
        ":lock",
    ],
    visibility = ["//visibility:public"],
)
