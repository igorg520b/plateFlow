#include "mainwindow.h"
#include <QApplication>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_sinks.h"

#include <QDebug>
#include <iostream>

#include <omp.h>


int main(int argc, char *argv[])
{
//    icy::CohesiveZone::CalculateAndPrintBMatrix();
    auto sharedFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("log.txt",true);
    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(std::make_shared<spdlog::sinks::stdout_sink_st>());
    sinks.push_back(sharedFileSink);
    auto combined_logger = std::make_shared<spdlog::logger>("name", begin(sinks), end(sinks));
    spdlog::register_logger(combined_logger);
    spdlog::set_default_logger(combined_logger);
    spdlog::flush_on(spdlog::level::info);

    spdlog::info("testing threads {}", omp_get_max_threads());
#pragma omp parallel
    {     spdlog::info("{}", omp_get_thread_num()); }
    std::cout << std::endl;

    spdlog::set_pattern("%v");

    QApplication a(argc, argv);
    QApplication::setApplicationName("plateFlow");
    QApplication::setApplicationVersion("1.1");

//    QSurfaceFormat fmt = QVTKOpenGLNativeWidget::defaultFormat();
//    QSurfaceFormat::setDefaultFormat(fmt);

    MainWindow w;
    w.resize(1400,900);
    w.show();
    return a.exec();
}
