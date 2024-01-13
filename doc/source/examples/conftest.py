def pytest_collectstart(collector):
    if collector.fspath and collector.fspath.ext == ".ipynb":
        collector.skip_compare += 'text/html', 'application/javascript', 'stderr',
