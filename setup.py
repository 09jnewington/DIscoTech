from setuptools import find_packages, setup

if __name__ == "__main__":
    setup(
        name="discotech",
        version="0.0.1",
        packages=find_packages(),
        author="Josh Newington",
        email="jn438@cam.ac.uk",
        description="A computational platform for natural project discovery and optimization.",
        package_data={
            "discotech": [
                "data/*",
            ],
        },
        include_package_data=True,
        entry_points={"console_scripts": ["discotech=discotech.__main__:main"]},
    )
