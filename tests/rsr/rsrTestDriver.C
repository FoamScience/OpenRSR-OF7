#define CATCH_CONFIG_RUNNER
#include "catch.H"
#include "error.H"

int main(int argc, char* argv[]) {

    // Cause FatalErrors and FatalIOErrors To Throw Exceptions
    Foam::FatalError.throwExceptions();
    Foam::FatalIOError.throwExceptions();

    // Run tests
    int result = Catch::Session().run(argc, argv);
    return result;
}
